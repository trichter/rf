"""
Frequency and time domain deconvolution.
"""
import numpy as np
from numpy import max, pi
from obspy.signal.util import next_pow_2
from scipy.fftpack import fft, ifft
from scipy.signal import correlate
try:
    from toeplitz import sto_sl
except ImportError:
    import warnings
    msg = 'Toeplitz import error. Time domain deconvolution will not work.'
    warnings.warn(msg)


def deconvolve(stream, method='time',
               source_components='LZ', response_components=None,
               set_tw='P',
               winsrc=None, tshift=None,
               **kwargs):
    """
    Deconvolve one component of a stream from all components.

    The deconvolutions are written to the data arrays of the stream. To keep
    the original data use the copy method of Stream.
    The stats dictionaries of the traces inside stream must have an 'onset'
    entry with a `~obspy.core.utcdatetime.UTCDateTime` object.
    This will be used for determining the data windows.

    :param stream: Stream object including response and source
    :param source_component: Name of component using for source function
    :param method:
        'time' -> use time domain deconvolution in `deconvt()`,\n
        'freq' -> use freqeuency domain deconvolution in `deconvf()`
    :type winsrc: tuple (start, end, taper)
    :param winsrc:
        data window for source function, in seconds relative to onset,\n
        default:\n
        (-10, 30, 5) for method='time',\n
        (-20, 80, 5) for method='freq'
    :type winrsp: tuple (start, end)
    :param winrsp: data window for response functions,\n
        just for method='time', default: (-20, 80)
    :type winrf: tuple (start, end)
    :param winrf: data window for results/deconvolution/receiver functions,\n
        just for method='time', default: (-20, 80)

    Other optional parameters are passed to the underlying deconvolution
    functions `deconvt()` and `deconvf()`
    .
    """
    if method not in ('time', 'freq'):
        raise NotImplementedError()
    # identify source and response components
    src = [tr for tr in stream if tr.stats.channel[-1] in source_components]
    if len(src) != 1:
        msg = 'Invalid number of source components. %d not equal to one.'
        raise ValueError(msg % len(src))
    src = src[0]
    rsp = [tr for tr in stream if response_components is None or
           tr.stats.channel[-1] in response_components]
    if not 0 < len(rsp) < 4:
        msg = 'Invalid number of response components. %d not between 0 and 4.'
        raise ValueError(msg % len(rsp))
    # define default time windows
    lenrsp = src.stats.endtime - src.stats.starttime
    onset = src.stats.onset - src.stats.starttime
    if set_tw == 'P' and method == 'time':
        winsrc = winsrc or (-10, 30, 5)
#        winrsp = winrsp or (-onset, lensec - onset)
#        winrf = winrf or (-onset, lensec - onset)
    elif set_tw == 'S' and method == 'time':
        winsrc = winsrc or (-10, 30, 5)
#        winrsp = winrsp or (onset - lensec, onset)
#        winrf = winrf or (onset - lensec, onset)
    elif set_tw == 'P':
        winsrc = winsrc or (-onset + 5, lenrsp - onset - 5, 5)
    elif set_tw == 'S':
        winsrc = winsrc or (onset - lenrsp + 5, onset - 5),

    winsrc = list(winsrc)
    if winsrc[0] < -onset:
        winsrc[0] = -onset
    if winsrc[1] - winsrc[0] > lenrsp:
        winsrc[1] = lenrsp - onset
    if tshift is None:
        tshift = -winsrc[0]

    # prepare source and response data and deconvolve
    sr = src.stats.sampling_rate
    onset_utc = src.stats.onset
    if src in rsp:
        src = src.copy()
    src.trim(onset_utc + winsrc[0], onset_utc + winsrc[1])
    src.taper(max_percentage=None, max_length=winsrc[2])
    rsp_data = [tr.data for tr in rsp]
    if method == 'time':
        lensrc = winsrc[1] - winsrc[0]
#        time_rf = winrf[1] - winrf[0]
#        shift = int(round(sr *
#                          ((winrsp[1] - winrsp[0] - winsrc[1] + winsrc[0] -
#                            time_rf) / 2 + winrsp[0] - winsrc[0] - winrf[0])))
#        length = int(round(time_rf * sr)) + 1
#        tshift = -winrf[0]
        shift = int(round(sr * (tshift - lensrc / 2)))
        length = int(round(lenrsp * sr)) + 1
        rf_data = deconvt(rsp_data, src.data, shift, length=length, **kwargs)
    else:
        rf_data = deconvf(rsp_data, src.data, sr, tshift=tshift, **kwargs)
    for tr in rsp:
        tr.stats.onset = tr.stats.onset + tshift + winsrc[0]
        tr.stats.tshift = tshift
    for i, tr in enumerate(rsp):
        tr.data = rf_data[i].real
    return rsp


def deconvf(rsp_list, src, sampling_rate, water=0.05, gauss=2., tshift=10.,
            pad=0, length=None, normalize=True, normalize_to_src=False,
            return_dict=False):
    """
    Frequency-domain deconvolution using waterlevel method.

    Deconvolve src from arrays in rsp_list.

    :param rsp_list: either a list of arrays containing the response functions
        or a single array
    :param src: array of source function
    :param sampling_rate: sampling rate of the data
    :param water: waterlevel to stabilize the deconvolution
    :param gauss: Gauss parameter of Low-pass filter
    :param tshift: delay time 0s will be at time tshift afterwards
    :param pad: multiply number of samples used for fft by 2**pad
    :param length: number of data points in results, optional
    :param normalize: if results are normalized
    :param normalize_to_src: True ->  normalized so that the maximum of a
        deconvolution of the source with itself is 1\n
        False -> normalized so that the maximum of the deconvolution of the
        first response array in rsp_list is 1
    :param return_dict: return additionally a lot of different parameters in a
        dict for debugging purposes

    :return: (list of) array(s) with deconvolution(s)
    """
    if length is None:
        length = len(src)
    N = length
    nfft = next_pow_2(N) * 2 ** pad
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
    if return_dict:
        if not normalize_to_src:
            spec_src = gauss * spec_src * spec_src_conj / spec_src_water
            rf_src = ifft(spec_src, nfft)[:N]
            norm = 1 / max(rf_src)
            rf_src = norm * rf_src
        ret_dict = {'rf_src': rf_src, 'rf_src_conj': spec_src_conj,
                    'spec_src_water': spec_src_water, 'freq': freq,
                    'gauss': gauss, 'norm': norm, 'N': N, 'nfft': nfft}
        return rf_list, ret_dict
    elif flag:
        return rf
    else:
        return rf_list


def _add_zeros(a, num, side='both'):
    """Add num zeros at side of array a"""
    return np.hstack([np.zeros(num)] * (side in ('both', 'left')) + [a] +
                     [np.zeros(num)] * (side in ('both', 'right')))


def _acorrt(a, num):
    """
    Not normalized auto-correlation of signal a.

    Sample 0 corresponds to zero lag time. Auto-correlation will consist of
    num samples. Correlation is performed in time domain by scipy.

    :param a: Data
    :param num: Number of returned data points
    :return: autocorrelation
    """
    return correlate(_add_zeros(a, num - 1, 'right'), a, 'valid')


def _xcorrt(a, b, num, zero_sample=0):
    """
    Not normalized cross-correlation of signals a and b.

    :param a,b: data
    :param num: The cross-correlation will consist of 2*num+1 samples.\n
        The sample with 0 lag time will be in the middle.
    :param zero_sample: Signals a and b are aligned around the middle of their
        signals.\n
        If zero_sample != 0 a will be shifted additionally to the left.
    :return: cross-correlation
    """
    if zero_sample != 0:
        a = _add_zeros(a, 2 * abs(zero_sample),
                       'left' if zero_sample > 0 else 'right')
    dif = len(a) - len(b) - 2 * num
    if dif > 0:
        b = _add_zeros(b, dif // 2)
    else:
        a = _add_zeros(a, -dif // 2)
    return correlate(a, b, 'valid')


def _toeplitz_real_sym(a, b):
    """
    Solve linear system Ax=b for real symmetric Toeplitz matrix A.

    :param a: first row of Toeplitz matrix A
    :param b: vector b
    :return: x=A^-1*b
    """
    return sto_sl(np.hstack((a, a[1:])), b, job=0)


# Gives similar results as a deconvolution with Seismic handler,
# but SH is faster
def deconvt(rsp_list, src, shift, spiking=1., length=None, normalize=True):
    """
    Time domain deconvolution.

    Deconvolve src from arrays in rsp_list.
    Calculate Toeplitz auto-correlation matrix of source, invert it, add noise
    and multiply it with cross-correlation vector of response and source.

    In one formula::

        RF = (STS + spiking*I)^-1 * STR

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
        STS = S^T*S = symmetric Toeplitz autocorrelation matrix
        STR = S^T*R = cross-correlation vector
        I... Identity

    :param rsp_list: either a list of arrays containing the response functions
        or a single array
    :param src: array of source function
    :param shift: shift the source by that amount of samples to the left side
        to get onset in RF at the desired time (negative -> shift source to the
        right side)\n
        shift = (middle of rsp window - middle of src window) +
        (0 - middle rf window)
    :param spiking: random noise added to autocorrelation (eg. 1.0, 0.1)
    :param length: number of data points in results
    :param normalize: normalize all results so that the maximum of the first
        result is 1

    :return: (list of) array(s) with deconvolution(s)
    """
    if length is None:
        length = len(src)
    flag = False
    RF_list = []
    STS = _acorrt(src, length)
    STS = STS / STS[0]
    STS[0] += spiking
    if not isinstance(rsp_list, (list, tuple)):
        flag = True
        rsp_list = [rsp_list]
    for rsp in rsp_list:
        STR = _xcorrt(rsp, src, length // 2, shift)
        if len(STR) > len(STS):
            STR = np.delete(STR, -1)
        RF = _toeplitz_real_sym(STS, STR)
        RF_list.append(RF)
    if normalize:
        norm = 1 / np.max(np.abs(RF_list[0]))
        for RF in RF_list:
            RF *= norm
    if flag:
        return RF
    else:
        return RF_list
