# Copyright 2013-2019 Tom Eulenfeld, MIT license
"""
Harmonic decomposition
"""
import numpy as np
import warnings
from rf.util import _add_processing_info


@_add_processing_info
def harmonics(stream, components='R', azim=0, method='time', **kwargs):
    """
    Perform harmonic decomposition of PRFs.

    This implements the method described in Bianchi et al 2010,
    doi:10.1029/2009JB007061 (equations in supplemental material) and in
    Park and Levin 2016, doi:10.1093/gji/ggw323 (equations 44/45 and 47).

    .. The equations in the two papers are equivalent, just written
       out slightly differently. Park and Levin have some sign errors.

    The harmonic components are returned in an array with rows for
    constant, cos, sin, cos2, and sin2 terms.
    The stats dictionaries of the traces in the input stream must have a
    'back_azimuth' entry and an 'event_time' entry.

    .. 'event_time' is used for trace sorting, under the assumption that
       it's a unique key and can be used to make sure that different
       components (e.g. R and T) are indexed in a common order.

    :param stream: RFStream including components to be decomposed
    :param components: names of components to use in decomposition;
        can be R, RT, Q, QT, or T (for PRFs)
    :param azim: azimuth along which to decompose the RFs (default: 0)
    :param method: domain for harmonic decomposition. Options are
        'time' -> time domain (Bianchi et al. 2010)
        'freq' -> frequency domain (Park and Levin 2016)
    :param \*\*kwargs: other kwargs are passed to underlying functions for
        building jacobians
    :return: harmonic components in a Stream

    .. note::
        The number of traces returned depends on the number of components.
        Two-component decomposition returns 10 harmonics: 5 modeled, 5
        unmodeled. One-component decomposition returns only 5 modeled harmonics
        for that single component.
    """
    # various checks for components and inputs
    if components not in ['RT', 'R', 'Q', 'QT', 'T']:  # TODO L/Z for SRFs?
        raise NotImplementedError('Component choice %s not supported' % components)
    for c in components:
        nrf = len(stream.select(channel='*%s' % c))
        if nrf == 0:
            raise ValueError('Component %s not in stream' % c)
        if nrf < 4:  # we probably shouldn't try to fit sinusoids with very few points
            warnings.warn('Not enough %s RFs for robust azimuthal fit' % c)
    if method not in ['time', 'freq']:
        raise NotImplementedError('Supported methods are time and freq')

    # short names to save some typing
    R = components[0]
    T = False
    if len(components) == 2:
        T = True        # and the second is T

    # get back azimuths (in degrees) with evt time as unique key
    baz = np.array([tr.stats.back_azimuth for tr in
                    stream.select(channel='*%s' % R).sort(['event_time'])])
    # get the RF trace data
    rfR = np.array([tr.data for tr in stream.select(channel='*%s' % R).sort(['event_time'])])
    if T:
        rfT = np.array([tr.data for tr in stream.select(channel='*T').sort(['event_time'])])

    # transform if working in frequency domain
    if method == 'freq':
        rfR = np.array([np.fft.fft(rfR[i]) for i in range(len(rfR))])
        if T:
            rfT = np.array([np.fft.fft(rfT[i]) for i in range(len(rfT))])

    # use methods to get the jacobian
    if T:
        jac = decomp_two(baz, azim=azim, **kwargs)
    else:
        jac = decomp_one(baz, azim=azim, **kwargs)

    # set up array for outputs
    npts = stream[0].stats.npts
    nharm = int(len(components)*5)
    if method == 'time':
        hd = np.zeros((nharm, npts))
    if method == 'freq':
        hd = np.zeros((nharm, npts), dtype=complex)

    # loop time points and fit at each one
    for i in range(npts):
        if not T:
            rfv = rfR[:, i]  # make data vector
        if T:
            rfv = np.hstack((rfR[:, i], rfT[:, i]))
        hd[:, i], _, _, _ = np.linalg.lstsq(jac, rfv, rcond=None)

    if method == 'freq':  # transform back if needed
        for i in range(nharm):
            hd[i] = np.fft.ifft(hd[i])
        hd = np.real(hd)

    # make the stream to return, with enough header info to get times later on
    out = stream[:1].copy()  # start a new stream with one trace
    for k in ['processing', 'event_id', 'event_depth', 'event_latitude', 'event_longitude',
              'event_magnitude', 'inclination', 'slowness', 'event_time', 'back_azimuth']:
        try:  # try to clean out some stats that no longer apply
            _ = out[0].stats.pop(k)
        except KeyError:
            pass
    out[0].stats['type'] = 'harmonic'
    out[0].data = hd[0]  # place harmonics in traces
    for i in range(1, nharm):
        new_trace = out[0].copy()
        new_trace.data = hd[i]
        out.append(new_trace)
    terms = ['constant', 'cos', 'sin', 'cos2', 'sin2']
    for i in range(5):  # add some metadata, replace channel and loc for sorting later
        out[i].stats['channel'] = str(i)
        out[i].stats['term'] = terms[i]
        out[i].stats['location'] = 'mod'
        out[i].stats['components'] = components
    if T:
        for i in range(5, 10):
            out[i].stats['channel'] = str(i-5)
            out[i].stats['term'] = terms[i-5]
            out[i].stats['location'] = 'unmod'
        out[i].stats['components'] = components
    return out


def decomp_two(baz, azim=0, scalars=(3, 0.3)):
    """ Build jacobian for harmonic decomposition with two components of RFs.

    :param baz: array of back azimuths for RFs in degrees
    :param azim: azimuth along which to decompose the RFs.
        This can be used with some kind of optimization minimize components
    :param scalars: scalar multipliers for R/Q and T constant terms.
        Park and Levin (2016) recommend 3 and 0.3 to bias toward radial component
        conversions, especially for datasets with patchy back-azimuthal coverage.
        (1,1) would weight Q/R and T equally in the regression.
    :return: jacobian matrix as np.ndarray, 10xN where N is the time series length
    """
    jacr = np.array([scalars[0]*np.ones(len(baz)),
                     np.cos(np.radians(baz-azim)),
                     np.sin(np.radians(baz-azim)),
                     np.cos(2*np.radians(baz-azim)),
                     np.sin(2*np.radians(baz-azim)),
                     np.zeros(len(baz)),
                     np.cos(np.radians(baz-azim)),
                     np.sin(np.radians(baz-azim)),
                     np.cos(2*np.radians(baz-azim)),
                     np.sin(2*np.radians(baz-azim))])
    jact = np.array([np.zeros(len(baz)),
                     -np.sin(np.radians(baz-azim)),
                     np.cos(np.radians(baz-azim)),
                     -np.sin(2*np.radians(baz-azim)),
                     np.cos(2*np.radians(baz-azim)),
                     scalars[1]*np.ones(len(baz)),
                     np.sin(np.radians(baz-azim)),
                     -np.cos(np.radians(baz-azim)),
                     np.sin(2*np.radians(baz-azim)),
                     -np.cos(2*np.radians(baz-azim))])
    jac = np.hstack((jacr, jact))
    return jac.T


def decomp_one(baz, azim=0):
    """ Build jacobian for harmonic decomposition for one RF component.

    :param baz: array of back azimuths for RFs in degrees
    :param azim: azimuth along which to decompose the RFs.
        This can be used with some kind of optimization minimize components
    :return: jacobian matrix as np.ndarray, 5xN where N is time series length
    """
    jac = np.array([np.ones(len(baz)),
                    np.cos(np.radians(baz-azim)),
                    np.sin(np.radians(baz-azim)),
                    np.cos(2*np.radians(baz-azim)),
                    np.sin(2*np.radians(baz-azim))])
    return jac.T
