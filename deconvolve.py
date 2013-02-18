
from numpy import max, pi
from obspy.signal.util import nextpow2
from scipy.fftpack import fft, ifft
from scipy.signal import correlate
import _toeplitz
import numpy as np

def deconv(stream, src_comp, method='time', **kwargs):
    if method == 'time':
        winsrc = kwargs.get('winsrc', (-10, 30, 5))
        winrsp = kwargs.get('winrsp', (-20, 80))
        winrf = kwargs.get('winrf', (-20, 80))
    elif method == 'freq':
        winsrc = kwargs.get('winsrc', (-20, 80, 5))
        tshift = kwargs.get('tshift', 10)
    else:
        raise Exception('Method is not a valid method')
    st = stream
    samp = st[0].stats.sampling_rate
    onset = st[0].stats.onset
    src = st.select(component=src_comp)[0]
    src_index = st.traces.index(src)
    src = src.copy()
    src.trim(onset + winsrc[0], onset + winsrc[1])
    src.taper(p=2 * winsrc[2] / (winsrc[1] - winsrc[0]))
    src = src.data
    rsp = [st[(i + src_index) % 3].data for i in range(3)]
    if method == 'time':
        time_rf = winrf[1] - winrf[0]
        shift = int(samp * ((winrsp[1] - winrsp[0] - winsrc[1] + winsrc[0] -
                         time_rf) / 2 + winrsp[0] - winsrc[0] - winrf[0]))
        length = int(time_rf * samp)
        rf_resp = deconvt(rsp, src, shift, length=length, **kwargs)
        tshift = -winrf[0]
    else:
        rf_resp = deconvf(rsp, src, samp, **kwargs)
    for i in range(3):
        st[(i + src_index) % 3].data = rf_resp[i].real
    for tr in st:
        tr.stats.tshift = tshift

def deconvf(rsp_list, src, sampling_rate, water=0.05, gauss=2., tshift=10.,
            pad=0, length=None, normalize=True, normalize_to_src=False,
            returnall=False):
    """
    Frequency-domain deconvolution using waterlevel method.

    rsp, src    data containing the response and
                source functions, respectively
    water       waterlevel to stabilize the deconvolution
    gauss       Gauss parameter of Low-pass filter
    tshift      shift the resulting function by that amount
    pad         multiply number of samples used for fft by 2**pad
    """
    if length == None:
        length = len(src)
    N = length
    nfft = nextpow2(N) * 2 ** pad
    freq = np.fft.fftfreq(nfft, d=1. / sampling_rate)
    gauss = np.exp(np.maximum(-(0.5 * 2 * pi * freq / gauss) ** 2, -700.) -
                   1j * tshift * 2 * pi * freq)

    spec_src = fft(src, nfft)
    spec_src_conj = np.conjugate(spec_src)
    spec_src_water = np.abs(spec_src * spec_src_conj)
    spec_src_water = np.maximum(spec_src_water, max(spec_src_water) * water)

    if normalize_to_src:
        spec_src = gauss * spec_src * spec_src_conj / spec_src_water
        rf_src = ifft(spec_src, nfft)[:N]
        #i1 = int((tshift-1)*sampling_rate)
        #i2 = int((tshift+1)*sampling_rate)
        norm = 1 / max(rf_src)
        rf_src = norm * rf_src

    flag = False
    if not isinstance(rsp_list, (list, tuple)):
        flag = True
        rsp_list = [rsp_list]
    rf_list = [ifft(gauss * fft(rsp, nfft) * spec_src_conj / spec_src_water,
                    nfft)[:N] for rsp in rsp_list]
    if normalize:
        if not normalize_to_src:
            norm = 1. / max(rf_list[0])
        for rf in rf_list:
            rf *= norm
    if returnall:
        if not normalize_to_src:
            spec_src = gauss * spec_src * spec_src_conj / spec_src_water
            rf_src = ifft(spec_src, nfft)[:N]
            norm = 1 / max(rf_src)
            rf_src = norm * rf_src
        return rf_list, rf_src, spec_src_conj, spec_src_water, freq, gauss, norm, N, nfft
    elif flag:
        return rf
    else:
        return rf_list

def add_zeros(a, num, side='both'):
    return np.hstack([np.zeros(num)] * (side in ('both', 'left')) + [a] +
                     [np.zeros(num)] * (side in ('both', 'right')))

def acorrt(a, num):
    """
    Return not normalized auto-correlation of signal a.
    
    Sample 0 corresponds to zero lag time. Auto-correlation will consist of
    num samples.
    """
    return correlate(add_zeros(a, num, 'right'), a, 'valid')

def xcorrt(a, b, num, zero_sample=0):
    """
    Return not normalized cross-correlation of signals a and b.
    
    zero_sample: Signals a and b are aligned around the middle of there signals.
                 If zero_sample != 0 a will be shifted additionally to the left.
    num: The cross-correlation will consist of 2*num+1 samples. The sample with
         0 lag time will be in the middle.
    """
    if zero_sample != 0:
        a = add_zeros(a, 2 * abs(zero_sample),
                      'left' if zero_sample > 0 else 'right')
    dif = len(a) - len(b) - 2 * num
    if dif > 0:
        b = add_zeros(b, dif // 2)
    else:
        a = add_zeros(a, -dif // 2)
    return correlate(a, b, 'valid')

def toeplitz(a1, b, a2=None):
    """
    Calculate inverse A^-1*b of Toeplitz matrix A.
    
    a1   first row of Toeplitz matrix A
    b    vector for multiplication with Inverse of Toeplitz matrix 
    """
    if len(a1.shape) == 1:
        a1 = a1[np.newaxis, :]
    if a2 == None:
        a2 = a1[:, 1:]
    elif len(a2.shape) == 1:
        a2 = a2[np.newaxis, :]
    if len(b.shape) == 1:
        bnew = b[np.newaxis, :]
    ret = _toeplitz.cbto_sl(a1, a2, bnew)
    if not a1.dtype == np.complex and not a2.dtype == np.complex:
        ret = np.real(ret)
    if len(b.shape) == 1:
        ret = ret[0, :]
    return ret

# Gives similar results as a deconvolution with Seismic handler,
# but SH is faster
def deconvt(rsp_list, src, shift, spiking=1, length=None, normalize=True):
    """
    Time domain deconvolution.

    Calculate Toeplitz auto-correlation matrix of source, invert it, add noise
    and multiply it with cross-correlation vector of response and source.

    In a formula:
    RF = (STS + spiking*I)^-1 * STR
    This function calculates RF.

    N... length
        ( S0   S-1  S-2 ... S-N+1 )
        ( S1   S0   S-1 ... S-N+2 )
    S = ( S2   ...                )
        ( ...                     )
        ( SN-1 ...          S0    )
    R = (R0 R1 ... RN-1)^T
    RF = (RF0 RF1 ... RFN-1)^T
    S... source matrix (shape N*N)
    R... response vector (length N)
    RF... receiver function (deconvolution) vector (length N)
    STS = S^T*S = Toeplitz autocorrelation matrix
    STR = S^T*R = cross-correlation vector
    I... Identity


    :parameters:
    rsp_list    (list of) data containing the response
    src         data of source
    spiking     random noise added to autocorrelation (eg. 1.0, 0.1)
    shift       shift the source by that amount of samples to the left side to get onset in RF at the right time
                (negative -> shift source to the right side)
                shift = (middle of rsp window - middle of src window) + (0 - middle rf window)
    length      number of data points of deconvolution
    normalize   if True normalize all deconvolutions so that the maximum of the
                first deconvolution is 1
    :return:    (list of) deconvolutions (length N)
    """
    if length == None:
        length = len(src)
    flag = False
    RF_list = []
    #STS = sito.xcorr.acorrt(src, length, demean=False, clipdata=False)
    STS = acorrt(src, length)
    STS = STS / STS[0]
    STS[0] += spiking
    #print shift, len(src), len(rsp_list), spiking, length
    if not isinstance(rsp_list, (list, tuple)):
        flag = True
        rsp_list = [rsp_list]
    for rsp in rsp_list:
        #STR = sito.xcorr.xcorrt(rsp, src, length // 2, shift, demean=False)
        STR = xcorrt(rsp, src, length // 2, shift)
        if len(STR) > len(STS):
            STR = np.delete(STR, -1)
        RF = toeplitz(STS, STR)
        RF_list.append(RF)
    if normalize:
        norm = 1 / np.max(np.abs(RF_list[0]))
        for RF in RF_list:
            RF *= norm
    if flag:
        return RF
    else:
        return RF_list

#def deconvf(stream, src_comp, water=0.01, gauss=2, tshift=10,
#            pad=0, winsrc=(-10, 30, 5), normalize=True):
#    """
#    Apply deconvolution in frequency domain.
#
#    :param water: waterlevel to stabilize the deconvolution (relative to data maximum)
#    :param gauss: Gauss parameter of averaging function (std of LP-filter)
#    :param tshift: shift the resulting function by that amount
#    :param pad: multiply number of samples used for fft by 2**pad.
#    :param window, start, end, where, lenslope: use only window (start,end) around
#        where of type window (with lenslope seconds of smoothing) of source function
#    :param return_real: just use the real part
#    """
#    st = stream
#    samp = st[0].stats.sampling_rate
#    onset = st[0].stats.onset
#    src = st.select(component=src_comp)[0]
#    src_index = st.index(src)
#    src = src.copy()
#    src.trim(onset + winsrc[0], onset + winsrc[1])
#    src.taper(p=2 * winsrc[2] / (winsrc[1] - winsrc[0]))
#    src = src.data
#    rsp = [st[(i + src_index) % 3].data for i in range(3)]
#    rf_resp = _deconvf(rsp, src, samp, water, gauss,
#                       tshift, pad, normalize=normalize)
#    for i in range(3):
#        st[(i + src_index) % 3].data = rf_resp[i].real
#    for tr in st:
#        tr.stats.tshift = tshift
#
#def deconvt(stream, src_comp, winsrc=(-20, 80, 5),
#            winrsp=(-20, 80), winrf=(-20, 80),
#            spiking=1, normalize=True):
#    """
#    Aply deconvolution in time-domain.
#    """
#    st = stream
#    samp = st[0].stats.sampling_rate
#    onset = st[0].stats.onset
#    src = st.select(component=src_comp).copy()
#    src.trim(onset + winsrc[0], onset + winsrc[1])
#    src.taper(p=2 * winsrc[2] / (winsrc[1] - winsrc[0]))
#    src = src[0].data
#    rsp = st[0].slice(onset + winrsp[0], onset + winrsp[1])
#    rsp = [tr.data for tr in rsp]
#    time_rf = winrf[1] - winrf[0]
#    shift = int(samp * ((winrsp[1] - winrsp[0] - winsrc[1] + winsrc[0] -
#                         time_rf) / 2 + winrsp[0] - winsrc[0] - winrf[0]))
#    rf_resp = _deconvt(rsp, src, spiking, shift, length=int(time_rf * samp),
#                         normalize=normalize)
#    st[0].data = rf_resp[0].real
#    st[1].data = rf_resp[1].real
#    st[2].data = rf_resp[2].real
#    for tr in st:
#        tr.stats.tshift = -winrf[0]

