# Copyright 2013-2019 Tom Eulenfeld, MIT license
"""
Frequency and time domain deconvolution.
"""
import numpy as np
from numpy import max, pi
from scipy.fftpack import fft, ifft, next_fast_len
from scipy.signal import correlate
from rf.util import _add_processing_info
from obspy.signal.util import next_pow_2
from copy import copy


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
        'iter' -> use iterative time domain deconvolution in `deconvi()`\n
        'func' -> user defined function (func keyword)
    :param func: Custom deconvolution function with the following signature

             def custom_deconv(rsp: RFStream, src: RFTrace, tshift=10,
                               **other_kwargs_possible) -> RFStream:

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
        functions `deconvt()`, `deconvf()`, and `deconvi()`

    .. note::
        If parameter normalize is not present in kwargs and source component is
        not excluded from the results by response_components, results will be
        normalized such that the maximum of the deconvolution of the trimmed
        and tapered source from the untouched source is 1. If the source is
        excluded from the results, the normalization will performed against
        the first trace in results.
    """
    if method not in ('time', 'freq', 'iter', 'func'):
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
    elif winsrc == 'P' and method == 'iter':
        winsrc = (-onset_sec, lenrsp_sec - onset_sec, 5)
    elif winsrc == 'S' and method == 'iter':  # TODO: test this
        winsrc = (-onset_sec, lenrsp_sec - onset_sec, 5)
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
    tshift = -winsrc[0]
    if method == 'time':
        shift = int(round(tshift * sr - len(src) // 2))
        rsp_data = [tr.data for tr in rsp]
        rf_data = deconvt(rsp_data, src.data, shift,  **kwargs)
        for i, tr in enumerate(rsp):
            tr.data = rf_data[i].real
    elif method == 'freq':
        rsp_data = [tr.data for tr in rsp]
        rf_data = deconvf(rsp_data, src.data, sr, tshift=tshift, **kwargs)
        for i, tr in enumerate(rsp):
            tr.data = rf_data[i].real
    elif method == 'iter':
        rsp_data = [tr.data for tr in rsp]
        rf_data, nit = deconvi(rsp_data, src.data, tshift, src.stats.delta, **kwargs)
        for i, tr in enumerate(rsp):
            tr.data = rf_data[i].real
            tr.stats['iterations'] = nit[i]
    else:
        rsp = func(stream.__class__(rsp), src, tshift=tshift, **kwargs)
    return stream.__class__(rsp)


def __get_length(rsp_list):
    if isinstance(rsp_list, (list, tuple)):
        rsp_list = rsp_list[0]
    return len(rsp_list)


def deconvf(rsp_list, src, sampling_rate, waterlevel=0.05, gauss=0.5,
            tshift=10., length=None, normalize=0, nfft=None,
            return_info=False):
    """
    Frequency-domain deconvolution using waterlevel method.

    Deconvolve src from arrays in rsp_list.

    :param rsp_list: either a list of arrays containing the response functions
        or a single array
    :param src: array with source function
    :param sampling_rate: sampling rate of the data
    :param waterlevel: waterlevel to stabilize the deconvolution
    :param gauss: Gauss parameter (standard deviation) of the
        Gaussian Low-pass filter,
        corresponds to cut-off frequency in Hz for a response value of
        exp(0.5)=0.607.
    :param tshift: delay time 0s will be at time tshift afterwards
    :param nfft: explicitely set number of samples for the fft
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
    if nfft is None:
        nfft = next_fast_len(N)
    freq = np.fft.fftfreq(nfft, d=1. / sampling_rate)
    gauss = np.exp(np.maximum(-0.5 * (freq / gauss) ** 2, -700.) -
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


# Gives similar results as a deconvolution with Seismic handler,
# but SH is faster
def deconvt(rsp_list, src, shift, spiking=1., length=None, normalize=0,
            solve_toeplitz='toeplitz'):
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
    :param solve_toeplitz: Python function to use as a solver for Toeplitz
        system. Possible values are
        ``'toeplitz'`` for using the toeplitz package (default),
        ``'scipy'`` for using solve_toeplitz from the SciPy package, or
        a custom function.

    :return: (list of) array(s) with deconvolution(s)
    """
    if solve_toeplitz == 'toeplitz':
        try:
            from toeplitz import sto_sl
        except ImportError:
            from warnings import warn
            warn('Toeplitz import error. Fallback to slower Scipy version.')
            from scipy.linalg import solve_toeplitz
        else:
            def solve_toeplitz(a, b):
                """
                Solve linear system Ax=b for real symmetric Toeplitz matrix A.

                :param a: first row of Toeplitz matrix A
                :param b: vector b
                :return: x=A^-1*b
                """
                return sto_sl(np.hstack((a, a[1:])), b, job=0)
    elif solve_toeplitz == 'scipy':
        from scipy.linalg import solve_toeplitz
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
        RF = solve_toeplitz(STS, STR)
        RF_list.append(RF)
    if normalize is not None:
        norm = 1 / np.max(np.abs(RF_list[normalize]))
        for RF in RF_list:
            RF *= norm
    if flag:
        return RF
    else:
        return RF_list

def _gauss_filter(dt, nft,f0):
    """
    Gaussian filter with width f0

    :param dt: sample spacing in seconds
    :param nft: length of filter in points
    :param f0: width of gaussian filter in seconds
    :return: array with Gaussian filter
    """
    f = np.fft.fftfreq(nft,dt)
    w = 2*pi*f
    gauss = np.exp(-0.25*(w/f0)**2)/dt
    return gauss

def _gfilter(x, nft, gauss, dt):
    """
    Apply a (Gaussian) filter to a data array

    :param x: array of data to filter
    :param nft: number of points for fft
    :param gauss: filter to apply in frequency domain, from _gauss_filter() 
        or elsewhere
    :param dt: sample spacing in seconds
    :return: real part of filtered array
    """
    xf = fft(x, n=nft)
    xnew = ifft(xf*gauss*dt, n=nft)
    return xnew.real

def _fft_correlate(a, b, nft):
    """
    Correlate two arrays; basically equivalent to summing two cross-correlations 
        with each shifted toward one end by half the array length

    :param a, b: data arrays
    :param nft: number of points for fft
    :return: array of real part of correlation
    """
    x = ifft(fft(a, n=nft) * np.conj(fft(b, n=nft)), n=nft)
    return x.real

def _phase_shift(x, nft, dt, tshift):
    """
    Shift array to account for time before onset

    :param x: array to shift
    :param nft: number of points for fft
    :param dt: sample spacing in seconds
    :param tshift: time to shift by in seconds
    :return: scaled real part of shifted array
    """
    xf = fft(x, n=nft)
    ish = int(tshift/dt)
    p = 2*pi*np.arange(1, nft+1)*ish/nft
    xf = xf*np.vectorize(complex)(np.cos(p), -np.sin(p))
    x = ifft(xf, n=nft)/np.cos(2*pi*ish/nft)
    return x.real

def deconvi(rsp, src, tshift, dt, f0=3.0, itmax=400, minderr=0.001, normalize=0):
    """
    Iterative deconvolution.

    Deconvolve src from arrays in rsp.
    Iteratively construct a spike train based on least-squares minimzation of the
    difference between one component of an observed seismogram and a predicted 
    signal generated by convolving the spike train with an orthogonal component
    of the seismogram.

    Reference: Ligorría, J. P., & Ammon, C. J. (1999). Iterative Deconvolution
    and Receiver-Function Estimation. Bulletin of the Seismological Society of
    America, 89, 5.

    :param rsp: either a list of arrays containing the response functions
        or a single array
    :param src: array with source function
    :param tshift: time window length before onset in seconds
    :param dt: sampling interval of data in seconds
    :param f0: width of Gaussian filter in seconds
    :param itmax: limit on number of iterations/spikes to add
    :param minderr: stop iteration when the change in error from adding another
        spike drops below this threshold
    :param normalize: normalize all results so that the maximum of the trace
        with the supplied index is 1. Set normalize to None for no normalization.

    :return: (list of) array(s) with deconvolution(s)
    """

    nt = len(src)       # number of points actually in trace
    ncomp = len(rsp)    # number of components we're looping over here

    #nfft = next_pow_2(nt)               # closest power of 2 for fft'ing the arrays
    nfft = nt
    RF_out = np.zeros((ncomp,nt))       # spike trains that we're going to make
    nit = np.zeros(ncomp)               # number of iterations each component uses

    for c in range(ncomp):  # loop over the responses
        rms = np.zeros(itmax)    # to store rms
        p0 = np.zeros(nfft)      # and rf for this component iteration

        r0 = np.pad(rsp[c],(0,nfft-nt))  # zero-pad the source and response arrays
        s0 = np.pad(src,(0,nfft-nt))     # (only matters if nfft!=nt for the sake of 2**)

        gaussF = _gauss_filter(dt, nfft, f0)  # construct and apply gaussian filter
        r_flt = _gfilter(r0, nfft, gaussF, dt)
        s_flt = _gfilter(s0, nfft, gaussF, dt)

        sft = fft(s0, nfft)  # fourier transform of the source
        rem_flt = copy(r_flt)  # thing to subtract from as spikes are added to p

        powerR = np.sum(r_flt**2)  # power in the response for scaling

        it = 0
        sumsq_i = 1
        d_error = 100*powerR + minderr
        maxlag = 0.5*nfft

        while np.abs(d_error) > minderr and it < itmax:  # loop iterations, add spikes
            rs = _fft_correlate(rem_flt, s_flt, nfft)  # correlate (what's left of) the num & demon
            rs = rs/np.sum(s_flt**2)  # scale the correlation

            i1 = np.argmax(np.abs(rs[0:int(maxlag) - 1]))  # index for getting spike amplitude
            # note that ^abs there means negative spikes are allowed
            amp = rs[i1]/dt

            p0[i1] = p0[i1] + amp  # add the amplitude of the spike to our spike-train RF
            p_flt = _gfilter(p0, nfft, gaussF, dt)  # gaussian filter the spike
            p_flt = _gfilter(p_flt, nfft, sft, dt)  # convolve with fft of source

            rem_flt = r_flt - p_flt  # subtract spike estimate from source to see what's left to model
            sumsq = np.sum(rem_flt**2)/powerR
            rms[it] = sumsq   # save rms
            d_error = 100*(sumsq_i - sumsq)  # check change in error as a result of this iteration

            sumsq_i = sumsq     # update rms
            it = it + 1         # and add one to the iteration count

        # once we get out of the loop:
        p_flt = _gfilter(p0, nfft, gaussF, dt)
        p_flt = _phase_shift(p_flt, nfft, dt, tshift)
        RF_out[c,:] = p_flt[:nt]  # save the RF for output
        nit[c] = it

    if normalize is not None:
        norm = 1 / np.max(np.abs(RF_out[normalize]))
        for RF in RF_out:
            RF *= norm

    return RF_out, nit
