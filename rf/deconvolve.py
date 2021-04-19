# Copyright 2013-2019 Tom Eulenfeld, MIT license
"""
Frequency and time domain deconvolution.
"""
from copy import copy
import numpy as np
from numpy import max, pi
from scipy.fftpack import fft, ifft, next_fast_len
from scipy.signal import correlate, detrend
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
        'time' -> use time domain deconvolution, see `deconv_time()`,\n
        'waterlevel' -> use frequency domain deconvolution with water level, see `deconv_waterlevel()`\n
        'iterative' -> use iterative time domain deconvolution, see `deconv_iterative()`\n
        'multitaper' -> use frequency domain multitaper deconvolution, see `deconv_multitaper()`\n
        'func' -> user defined function (func keyword)
    :param func: Custom deconvolution function with the following signature

             def deconv_custom(rsp: RFStream, src: RFTrace, tshift=10,
                               \*\*other_kwargs_possible) -> RFStream:

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
        functions

    .. note::
        If parameter normalize is not present in kwargs and source component is
        not excluded from the results by response_components, results will be
        normalized such that the maximum of the deconvolution of the trimmed
        and tapered source from the untouched source is 1. If the source is
        excluded from the results, the normalization will performed against
        the first trace in results.

    .. note::
        If multitaper deconvolution is used and a stream of (pre-event) noise
        is not present in kwargs, noise will be sampled from the data in a
        pre-event window whose length depends on the trace length prior to the
        onset time.
    """
    if method == 'freq':
        method = 'waterlevel'
    if method not in ('time', 'waterlevel', 'iterative', 'multitaper', 'func'):
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
    if winsrc == 'P' and method in ('time', 'multitaper', 'func'):
        winsrc = (-10, 30, 5)
    elif winsrc == 'S' and method in ('time','multitaper', 'func'):  # TODO: test this
        winsrc = (-10, 30, 5)
    elif winsrc == 'P' and method == 'iterative':
        winsrc = (-onset_sec, lenrsp_sec - onset_sec, 0)
    elif winsrc == 'S' and method == 'iterative':  # TODO: test this
        winsrc = (-onset_sec, lenrsp_sec - onset_sec, 0)
    elif winsrc == 'P' and method == 'waterlevel':
        winsrc = (-onset_sec, lenrsp_sec - onset_sec, 5)
    elif winsrc == 'S' and method == 'waterlevel':
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
        rf_data = deconv_time(rsp_data, src.data, shift,  **kwargs)
        for i, tr in enumerate(rsp):
            tr.data = rf_data[i].real
    elif method == 'waterlevel':
        rsp_data = [tr.data for tr in rsp]
        rf_data = deconv_waterlevel(rsp_data, src.data, sr, tshift=tshift,
                                    **kwargs)
        for i, tr in enumerate(rsp):
            tr.data = rf_data[i].real
    elif method == 'iterative':
        rsp_data = [tr.data for tr in rsp]
        rf_data, nit, _ = deconv_iterative(rsp_data, src.data, sr,
                                           tshift=tshift, **kwargs)
        for i, tr in enumerate(rsp):
            tr.data = rf_data[i].real
            tr.stats['iterations'] = nit[i]
    elif method == 'multitaper':
        noise = kwargs.pop('noise',None)
        if noise is None:  # no kwarg, grab from pre-event time series
            onset_rsp = rsp[0].stats.onset - rsp[0].stats.starttime
            noise = stream.copy().trim2(-onset_rsp,-5,'onset') # NOTE window length will vary
        # noise is not None (kwarg provided), noise should be Stream() or RFStream()
        # so now we grab that data and make a list of arrays
        nse_data = [tr.data for tr in noise if response_components is None or
                    tr.stats.channel[-1] in response_components]
        rsp_data = [tr.data for tr in rsp]
        rf_data = deconv_multitaper(rsp_data, src.data, nse_data, sr, -tshift,
                                    **kwargs)
        for i, tr in enumerate(rsp):
            tr.data = rf_data[i].real
    else:
        rsp = func(stream.__class__(rsp), src, tshift=tshift, **kwargs)
    return stream.__class__(rsp)


def __get_length(rsp_list):
    if isinstance(rsp_list, (list, tuple)):
        rsp_list = rsp_list[0]
    return len(rsp_list)


def deconv_waterlevel(rsp_list, src, sampling_rate, waterlevel=0.05, gauss=0.5,
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
    dt = 1. / sampling_rate
    ffilt = _phase_shift_filter(nfft, dt, tshift)
    if gauss is not None:
        ffilt = _gauss_filter(dt, nfft, gauss, waterlevel=-700) * ffilt
    spec_src = fft(src, nfft)
    spec_src_conj = np.conjugate(spec_src)
    spec_src_water = np.abs(spec_src * spec_src_conj)
    spec_src_water = np.maximum(
        spec_src_water, max(spec_src_water) * waterlevel)

    if normalize == 'src':
        spec_src = ffilt * spec_src * spec_src_conj / spec_src_water
        rf_src = ifft(spec_src, nfft)[:N]
        norm = 1 / max(rf_src)
        rf_src = norm * rf_src

    flag = False
    if not isinstance(rsp_list, (list, tuple)):
        flag = True
        rsp_list = [rsp_list]
    rf_list = [ifft(ffilt * fft(rsp, nfft) * spec_src_conj / spec_src_water,
                    nfft)[:N] for rsp in rsp_list]
    if normalize not in (None, 'src'):
        norm = 1. / max(rf_list[normalize])
    if normalize is not None:
        for rf in rf_list:
            rf *= norm
    if return_info:
        if normalize not in (None, 'src'):
            spec_src = ffilt * spec_src * spec_src_conj / spec_src_water
            rf_src = ifft(spec_src, nfft)[:N]
            norm = 1 / max(rf_src)
            rf_src = norm * rf_src
        info = {'rf_src': rf_src, 'rf_src_conj': spec_src_conj,
                'spec_src_water': spec_src_water,
                'freq': np.fft.fftfreq(nfft, d=dt),
                'gauss': ffilt, 'norm': norm, 'N': N, 'nfft': nfft}
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
def deconv_time(rsp_list, src, shift, spiking=1., length=None, normalize=0,
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


def _gauss_filter(dt, nft, f0, waterlevel=None):
    """
    Gaussian filter with width f0

    :param dt: sample spacing in seconds
    :param nft: length of filter in points
    :param f0: Standard deviation of the Gaussian Low-pass filter,
        corresponds to cut-off frequency in Hz for a response value of
        exp(0.5)=0.607.
    :param waterlevel: waterlevel for eliminating very low values
        (default: no waterlevel)
    :return: array with Gaussian filter frequency response
    """
    f = np.fft.fftfreq(nft, dt)
    gauss_arg = -0.5 * (f/f0) ** 2
    if waterlevel is not None:
        gauss_arg = np.maximum(gauss_arg, waterlevel)
    return np.exp(gauss_arg)


def _phase_shift_filter(nft, dt, tshift):
    """
    Construct filter to shift an array to account for time before onset

    :param nft: number of points for fft
    :param dt: sample spacing in seconds
    :param tshift: time to shift by in seconds
    :return: shifted array
    """
    freq = np.fft.fftfreq(nft, d=dt)
    return np.exp(-2j * pi * freq * tshift)


def _apply_filter(x, filt):
    """
    Apply a filter defined in frequency domain to a data array

    :param x: array of data to filter
    :param filter: filter to apply in frequency domain,
        e.g. from _gauss_filter()
    :return: real part of filtered array
    """
    nfft = len(filt)
    xf = fft(x, n=nfft)
    xnew = ifft(xf*filt, n=nfft)
    return xnew.real


def deconv_iterative(rsp, src, sampling_rate, tshift=10, gauss=0.5, itmax=400,
                     minderr=0.001, mute_shift=False, normalize=0):
    """
    Iterative deconvolution.

    Deconvolve src from arrays in rsp.
    Iteratively construct a spike train based on least-squares minimzation of the
    difference between one component of an observed seismogram and a predicted
    signal generated by convolving the spike train with an orthogonal component
    of the seismogram.

    Reference: Ligorria, J. P., & Ammon, C. J. (1999). Iterative Deconvolution
    and Receiver-Function Estimation. Bulletin of the Seismological Society of
    America, 89, 5.

    :param rsp: either a list of arrays containing the response functions
        or a single array
    :param src: array with source function
    :param sampling_rate: sampling rate of the data
    :param tshift: delay time 0s will be at time tshift afterwards
    :param gauss: Gauss parameter (standard deviation) of the
        Gaussian Low-pass filter,
        corresponds to cut-off frequency in Hz for a response value of
        exp(0.5)=0.607.
    :param itmax: limit on number of iterations/spikes to add
    :param minderr: stop iteration when the change in error from adding another
        spike drops below this threshold
    :param mute_shift: Mutes all samples at beginning of trace
        (lenght given by time shift).
        For `len(src)==len(rsp)` this mutes all samples before the onset.
    :param normalize: normalize all results so that the maximum of the trace
        with the supplied index is 1. Set normalize to None for no normalization.

    :return: (list of) array(s) with deconvolution(s)
    """

    ncomp = len(rsp)    # number of components we're looping over here
    dt = 1 / sampling_rate
    RF_out = []
    it_out = []  # number of iterations each component uses
    rms_out = []

    for c in range(ncomp):  # loop over the responses
        r0 = rsp[c]
        nfft = next_fast_len(2 * len(r0))
        rms = np.zeros(itmax)    # to store rms
        p0 = np.zeros(nfft)      # and rf for this component iteration

        gaussF = _gauss_filter(dt, nfft, gauss)  # construct and apply gaussian filter
        r_flt = _apply_filter(r0, gaussF)
        s_spec_flt = fft(src, nfft) * gaussF  # spectrum of source
        powerS = np.sum(ifft(s_spec_flt).real ** 2)  # power in the source for scaling
        powerR = np.sum(r_flt**2)  # power in the response for scaling
        rem_flt = copy(r_flt)  # thing to subtract from as spikes are added to p0

        it = 0
        sumsq_i = 1
        d_error = 100*powerR + minderr
        mute_min = len(r0)
        mute_max = nfft if mute_shift else nfft - int(tshift*sampling_rate)
        while np.abs(d_error) > minderr and it < itmax:  # loop iterations, add spikes
            rs = ifft(fft(rem_flt) * np.conj(s_spec_flt)).real  # correlate (what's left of) the num & demon, scale
            rs = rs / powerS / dt  # scale the correlation
            rs[mute_min:mute_max] = 0
            i1 = np.argmax(np.abs(rs))  # index for getting spike amplitude
            # note that ^abs there means negative spikes are allowed
            p0[i1] = p0[i1] + rs[i1]  # add the amplitude of the spike to our spike-train RF
            p_flt = ifft(fft(p0) * s_spec_flt).real * dt  # convolve with source

            rem_flt = r_flt - p_flt  # subtract spike estimate from source to see what's left to model
            sumsq = np.sum(rem_flt**2) / powerR
            rms[it] = sumsq   # save rms
            d_error = 100 * (sumsq_i - sumsq)  # check change in error as a result of this iteration

            sumsq_i = sumsq     # update rms
            it = it + 1         # and add one to the iteration count
        # once we get out of the loop:
        shift_filt = _phase_shift_filter(nfft, dt, tshift)
        p_flt = _apply_filter(p0, gaussF * shift_filt)
        RF_out.append(p_flt[:len(r0)])  # save the RF for output
        it_out.append(it)
        rms_out.append(rms)

    if normalize is not None:
        norm = 1 / np.max(np.abs(RF_out[normalize]))
        for RF in RF_out:
            RF *= norm

    return RF_out, it_out, rms_out


def deconv_multitaper(rsp, src, nse, sampling_rate, tshift, gauss=0.5,
                      K=3, tband=4, T=10, olap=0.75, normalize=0):
    """
    Multitaper frequency domain deconvolution

    Deconvolve src from arrays in rsp.

    Based on mtdecon.f by G. Helffrich with some slight modifications

    References:
    Helffrich, G (2006). Extended-time multitaper frequency domain
    cross-correlation receiver function estimation. Bulletin of the
    Seismological Society of America, 96 (1).
    Shibutani, T., Ueno, T., & Hirahara, K. (2008). Improvement in the
    Extended-Time Multitaper Receiver Function Estimation Technique.
    Bulletin of the Seismological Society of America, 98 (2).

    :param rsp: a list of arrays containing the response functions
    :param src: array of source function
    :param nse: a list of arrays containing samples of pre-event noise
    :param sampling_rate: sampling rate of the data
    :param tshift: shift the source by that amount of samples to the left side
        to get onset in RF at the desired time (src time window length pre-onset)
    :param gauss: Gauss parameter (standard deviation) of the
        Gaussian Low-pass filter,
        corresponds to cut-off frequency in Hz for a response value of
        exp(0.5)=0.607.
    :param K: number of Slepian tapers to use (default: 3)
    :param tband: time-bandwidth product (default: 4)
    :param T: time length of taper window in seconds (default: 10)
    :param olap: window overlap between 0 and 1 (default: 0.75)
    :param normalize: normalize all results so that the maximum of the trace
        with supplied index is 1. Set normalize to None for no normalization.

    :return: (list of) array(s) with deconvolution(s)
    """
    try:
        import mtspec as mt
    except ImportError as ex:
        msg = 'mtspec package is needed for multitaper deconvolution'
        raise ImportError(msg) from ex

    # check src trace length < rsp trace length:
    assert len(src) < len(rsp[0]), 'source wavelet must be shorter than response'
    assert len(nse[0]) <= len(rsp[0]), 'noise should not be longer than response'

    nft = len(rsp[0])  # length of final arrays
    dt = 1./sampling_rate  # sample spacing

    # calculate multitapers with mtspec
    ntap = int(round(T/dt))  # #points in each taper

    tap, el, _ = mt.multitaper.dpss(ntap, tband, K); tap = tap.T
    if K >= 2:
        tap[1] = -tap[1]  # adjust to match sign convention
    if K >= 3:
        tap[2] = -tap[2]
    if K > 3:
        print('Warning: taper signs may need adjustment for K>3')

    ofac = 1 / (1 - olap)  # calculate overlap factor

    sfft = np.zeros((K,nft), dtype=np.complex)  # array for holding freq-domain source estimate
    src = detrend(src, type='constant')  # demean and detrend source
    src = detrend(src, type='linear')

    # window and taper the source
    nwav = len(src)   # number of nonzero points for the source
    nwin = int(ofac*nwav/ntap)  # number of windows needed to cover the source

    src_pad = np.zeros(nft)
    src_pad[:nwav] = src
    for j in range(nwin):
        js = int(j*ntap/ofac)  # get indices for overlapping taper window
        je = int(js + ntap)
        for k in range(K):
            bit = np.zeros(nft)
            bit[js:je] = tap[k] # insert window taper into full-length zeros
            bit = bit*src_pad   # window/taper the padded source trace
            bit = fft(bit)      # transform tapered bit to frequency domain
            sfft[k] = sfft[k] + bit  # sum to total freq-domain source estimate

    # multiply by factor to account for short wavelet relative to fft
    fac = nwav/float(nwin*nft)
    sfft = fac*sfft

    # loop response components and do a similar window/taper/transform thing
    ncomp = len(rsp)
    RF_out = np.zeros((ncomp, nft))
    for c in range(ncomp):
        dat = rsp[c]; nos = nse[c]  # pick out component
        dfft = np.zeros((K,nft), dtype=np.complex)  # arrays for storing freq-domain estimates
        nfft = np.zeros((K,nft), dtype=np.complex)

        # window, taper, and transform the noise
        nos = detrend(nos, type='constant')
        nos = detrend(nos, type='linear')
        nos_pad = np.zeros(nft)
        nos_pad[:len(nos)] = nos
        nwin = int(ofac*len(nos)/ntap)
        for j in range(nwin):
            js = int(j*ntap/ofac)
            je = int(js + ntap)
            if je < nft:
                for k in range(K):
                    bit = np.zeros(nft)
                    bit[js:je] = tap[k]
                    bit = bit*nos_pad
                    bit = fft(bit)
                    nfft[k] = nfft[k] + bit

        fac = nwav/float(nwin*nft)  # scale for noise trace length
        nfft = nfft*fac

        # calculate power in the noise window
        s0 = np.zeros(nft)
        for k in range(K):
            s0 = s0 + (np.real(nfft[k])**2 + np.imag(nfft[k])**2)/el[k]

        # window, taper, and transform the response traces
        dat = detrend(dat, type='constant')
        dat = detrend(dat, type='linear')
        nwin = int(ofac*nft/ntap)-1
        for j in range(nwin):
            js = int(j*ntap/ofac)
            je = int(js + ntap)
            if je < nft:
                for k in range(K):
                    bit = np.zeros(nft)
                    bit[js:je] = tap[k]
                    bit = bit*dat
                    bit = fft(bit)
                    dfft[k] = dfft[k] + bit

        fac = float(nft)/float(nwin*nft)
        dfft = dfft*fac

        # deconvolve
        num = np.zeros(nft)
        denom = np.zeros(nft)
        for k in range(K):
            num = num + sfft[k]*np.conj(dfft[k])
            denom = denom + sfft[k]*np.conj(sfft[k])
        recF = num/(denom + s0)

        # time-shift for onset and apply gaussian lowpass
        if gauss is not None:
            recF = recF * _gauss_filter(dt, nft, gauss)
        recF = recF * _phase_shift_filter(nft, dt, tshift-dt)  # one sample gets lost in the shuffle

        # inverse fft, put in output array
        recT = ifft(recF)
        RF_out[c] = copy(recT.real[::-1])

    if normalize is not None:
        norm = 1 / np.max(np.abs(RF_out[normalize]))
        for RF in RF_out:
            RF *= norm

    return RF_out
