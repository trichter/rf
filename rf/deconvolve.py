# Copyright 2013-2016 Tom Eulenfeld, MIT license
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

from rf.util import _add_processing_info


def __find_nearest(array, value):
    """http://stackoverflow.com/a/26026189"""
    idx = np.searchsorted(array, value, side='left')
    expr = np.abs(value - array[idx - 1]) < np.abs(value - array[idx])
    if idx > 0 and (idx == len(array) or expr):
        return idx - 1
    else:
        return idx


@_add_processing_info
def deconvolve(stream, method='time', func=None,
               source_components='LZ', response_components=None,
               winsrc='P', **kwargs):
    """
    Deconvolve one component of a stream from other components.

    The deconvolutions are written to the data arrays of the stream. To keep
    the original data use the copy method of Stream.
    The stats dictionaries of the traces inside stream must have an 'onset'
    entry with a `~obspy.core.utcdatetime.UTCDateTime` object.
    This will be used for determining the data windows.

    :param stream: stream including responses and source
    :param method:
        'time' -> use time domain deconvolution in `deconvt()`,\n
        'freq' -> use frequency domain deconvolution in `deconvf()`\n
        'func' -> user defined function (func keyword)
    :param source_components: names of components identifying the source traces,
        e.g. 'LZ' for P receiver functions and 'QR' for S receiver functions
    :param response_components: names of components identifying the response
        traces (default None: all traces are used as response)
    :type winsrc: tuple (start, end, taper)
    :param winsrc: data window for source function, in seconds relative to
        onset. The function will be cosine tapered at both ends (seconds).\n
        winsrc can also be a string ('P' or 'S'). In this case the function
        defines a source time window appropriate for this type of receiver
        function and deconvolution method (see source code for details).
    :param \*\*kwargs: other kwargs are passed to the underlying deconvolution
        functions `deconvt()` and `deconvf()`

    .. note::
        If parameter normalize is not present in kwargs and source component is
        not excluded from the results by response_components, results will be
        normalized such that the maximum of the deconvolution of the trimmed
        and tapered source from the untouched source is 1. If the source is
        excluded from the results, the normalization will performed against
        the first trace in results.
    """
    if method not in ('time', 'freq', 'func'):
        raise NotImplementedError()
    # identify source and response components
    src = [tr for tr in stream if tr.stats.channel[-1] in source_components]
    if len(src) != 1:
        msg = 'Invalid number of source components. %d not equal to one.'
        raise ValueError(msg % len(src))
    src = src[0]
    rsp = [tr for tr in stream if response_components is None or
           tr.stats.channel[-1] in response_components]
    if 'normalize' not in kwargs and src in rsp:
        kwargs['normalize'] = rsp.index(src)
    if not 0 < len(rsp) < 4:
        msg = 'Invalid number of response components. %d not between 0 and 4.'
        raise ValueError(msg % len(rsp))

    sr = src.stats.sampling_rate
    # shift onset to time of nearest data sample to circumvent complications
    # for data with low sampling rate and method='time'
    idx = __find_nearest(src.times(), src.stats.onset - src.stats.starttime)
    src.stats.onset = onset = src.stats.starttime + idx * src.stats.delta
    for tr in rsp:
        tr.stats.onset = onset
    # define default time windows
    lenrsp_sec = src.stats.endtime - src.stats.starttime
    onset_sec = onset - src.stats.starttime
    if winsrc == 'P' and method == 'time':
        winsrc = (-10, 30, 5)
    elif winsrc == 'S' and method == 'time':
        winsrc = (-10, 30, 5)
    elif winsrc == 'P':
        winsrc = (-onset_sec, lenrsp_sec - onset_sec, 5)
    elif winsrc == 'S':
        winsrc = (-10, lenrsp_sec - onset_sec, 5)
#    winsrc = list(winsrc)
#    if winsrc[0] < -onset_sec:
#        winsrc[0] = -onset_sec
#    if winsrc[1] > lenrsp_sec - onset_sec:
#        winsrc[1] = lenrsp_sec - onset_sec
    # prepare source and response list
    if src in rsp:
        src = src.copy()
    src.trim(onset + winsrc[0], onset + winsrc[1], pad=True, fill_value=0.)
    src.taper(max_percentage=None, max_length=winsrc[2])
    rsp_data = [tr.data for tr in rsp]
    tshift = -winsrc[0]
    if method == 'time':
        shift = int(round(tshift * sr - len(src) // 2))
        rf_data = deconvt(rsp_data, src.data, shift,  **kwargs)
    elif method == 'freq':
        rf_data = deconvf(rsp_data, src.data, sr, tshift=tshift, **kwargs)
    else:
        rf_data = func(rsp_data, src.data, sr=sr, tshift=tshift, **kwargs)
    for i, tr in enumerate(rsp):
        tr.data = rf_data[i].real
    return stream.__class__(rsp)


def __get_length(rsp_list):
    if isinstance(rsp_list, (list, tuple)):
        rsp_list = rsp_list[0]
    return len(rsp_list)


def deconvf(rsp_list, src, sampling_rate, waterlevel=0.05, gauss=2.,
            tshift=10., pad=0, length=None, normalize=0, return_info=False):
    """
    Frequency-domain deconvolution using waterlevel method.

    Deconvolve src from arrays in rsp_list.

    :param rsp_list: either a list of arrays containing the response functions
        or a single array
    :param src: array with source function
    :param sampling_rate: sampling rate of the data
    :param waterlevel: waterlevel to stabilize the deconvolution
    :param gauss: Gauss parameter of Low-pass filter
    :param tshift: delay time 0s will be at time tshift afterwards
    :param pad: multiply number of samples used for fft by 2**pad
    :param length: number of data points in results, optional
    :param normalize: normalize all results so that the maximum of the trace
        with supplied index is 1. Set normalize to 'src' to normalize
        for the maximum of the prepared source. Set normalize to None for no
        normalization.
    :param return_info: return additionally a lot of different parameters in a
        dict for debugging purposes

    :return: (list of) array(s) with deconvolution(s)
    """
    if length is None:
        length = __get_length(rsp_list)
    N = length
    nfft = next_pow_2(N) * 2 ** pad
    freq = np.fft.fftfreq(nfft, d=1. / sampling_rate)
    gauss = np.exp(np.maximum(-(0.5 * 2 * pi * freq / gauss) ** 2, -700.) -
                   1j * tshift * 2 * pi * freq)

    spec_src = fft(src, nfft)
    spec_src_conj = np.conjugate(spec_src)
    spec_src_water = np.abs(spec_src * spec_src_conj)
    spec_src_water = np.maximum(
        spec_src_water, max(spec_src_water) * waterlevel)

    if normalize == 'src':
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
    if normalize not in (None, 'src'):
        norm = 1. / max(rf_list[normalize])
    if normalize is not None:
        for rf in rf_list:
            rf *= norm
    if return_info:
        if normalize not in (None, 'src'):
            spec_src = gauss * spec_src * spec_src_conj / spec_src_water
            rf_src = ifft(spec_src, nfft)[:N]
            norm = 1 / max(rf_src)
            rf_src = norm * rf_src
        info = {'rf_src': rf_src, 'rf_src_conj': spec_src_conj,
                'spec_src_water': spec_src_water, 'freq': freq,
                'gauss': gauss, 'norm': norm, 'N': N, 'nfft': nfft}
        return rf_list, info
    elif flag:
        return rf
    else:
        return rf_list


def _add_zeros(a, numl, numr):
    """Add zeros at left and rigth side of array a"""
    return np.hstack([np.zeros(numl), a, np.zeros(numr)])


def _acorrt(a, num):
    """
    Not normalized auto-correlation of signal a.

    Sample 0 corresponds to zero lag time. Auto-correlation will consist of
    num samples. Correlation is performed in time domain with scipy.

    :param a: Data
    :param num: Number of returned data points
    :return: autocorrelation
    """
    return correlate(_add_zeros(a, 0, num - 1), a, 'valid')


def _xcorrt(a, b, num, zero_sample=0):
    """
    Not normalized cross-correlation of signals a and b.

    :param a,b: data
    :param num: The cross-correlation will consist of num samples.\n
        The sample with 0 lag time will be in the middle.
    :param zero_sample: Signals a and b are aligned around the middle of their
        signals.\n
        If zero_sample != 0 a will be shifted additionally to the left.
    :return: cross-correlation
    """
    if zero_sample > 0:
        a = _add_zeros(a, 2 * abs(zero_sample), 0)
    elif zero_sample < 0:
        a = _add_zeros(a, 0, 2 * abs(zero_sample))
    dif = len(a) - len(b) + 1 - num
    if dif > 0:
        b = _add_zeros(b, (dif + 1) // 2, dif // 2)
    else:
        a = _add_zeros(a, (-dif + 1) // 2, (-dif) // 2)
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
def deconvt(rsp_list, src, shift, spiking=1., length=None, normalize=0):
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
    :param normalize: normalize all results so that the maximum of the trace
        with supplied index is 1. Set normalize to None for no normalization.

    :return: (list of) array(s) with deconvolution(s)
    """
    if length is None:
        length = __get_length(rsp_list)
    flag = False
    RF_list = []
    STS = _acorrt(src, length)
    STS = STS / STS[0]
    STS[0] += spiking
    if not isinstance(rsp_list, (list, tuple)):
        flag = True
        rsp_list = [rsp_list]
    for rsp in rsp_list:
        STR = _xcorrt(rsp, src, length, shift)
        assert len(STR) == len(STS)
        RF = _toeplitz_real_sym(STS, STR)
        RF_list.append(RF)
    if normalize is not None:
        norm = 1 / np.max(np.abs(RF_list[normalize]))
        for RF in RF_list:
            RF *= norm
    if flag:
        return RF
    else:
        return RF_list
