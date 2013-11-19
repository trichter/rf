"""
Classes and functions for receiver function calculation.
"""

from math import pi, sin
from operator import attrgetter
import warnings

import numpy as np
from obspy.core import Stream, Trace
from obspy.core.util import AttribDict
from obspy.core.util.geodetics import gps2DistAzimuth, kilometer2degrees
from obspy.taup.taup import getTravelTimes
from rf.deconvolve import deconv
from rf.simple_model import load_model


def __get_event_attr(h1, h2):
    return lambda event: event[h1][0][h2]


def __rel2UTC(stats, head):
    return stats.starttime + stats[head]


def __UTC2rel(stats, head):
    return stats[head] - stats.starttime


STATION_GETTER = (('station_latitude', attrgetter('latitude')),
                  ('station_longitude', attrgetter('longitude')),
                  ('station_elevation', attrgetter('elevation')))
# TODO: map event_id
EVENT_GETTER = (('event_latitude', __get_event_attr('origins', 'latitude')),
                ('event_longitude', __get_event_attr('origins', 'longitude')),
                ('event_depth', __get_event_attr('origins', 'depth')),
                ('event_magnitude', __get_event_attr('magnitudes', 'mag')),
                ('event_time', __get_event_attr('origins', 'time')))
HEADERS = zip(*STATION_GETTER)[0] + zip(*EVENT_GETTER)[0] + (
    'onset', 'distance', 'back_azimuth', 'inclination', 'slowness')
FORMATHEADERS = {'sac': ('stla', 'stlo', 'stel', 'evla', 'evlo', 'evdp', 'mag',
                         'o', 'a', 'gcarc', 'baz', 'user0', 'user1'),
                 # fields 'dcvreg', 'dcvinci' and 'pwdw' are violated for
                 # station information
                 'sh': ('DCVREG', 'DCVINCI', 'PWDW', 'LAT', 'LON', 'DEPTH',
                        'MAGNITUDE', 'ORIGIN', 'P-ONSET', 'DISTANCE',
                        'AZIMUTH', 'INCI', 'SLOWNESS')}
_HEADER_CONVERSIONS = {'sac': {'onset': (__rel2UTC, __UTC2rel),
                               'event_time': (__rel2UTC, __UTC2rel)}}
_HEADERS_EXAMPLE = (50.3, -100.2, 400.3, -20.32, 10., 12.4, 6.5, -40.432,
                    20.643, 57.6, 90.1, 10.2, 10.)


class RFStream(Stream):
    """
    Class providing the necessary functions for receiver function calculation.

    To initialize a RFStream from a ObsPy stream use

    >>> rfstream = RFStream(stream=obspy_stream)

    Format specific headers are loaded into the stats object of all traces.
    """
    def __init__(self, traces=None, stream=None):
        if stream is not None:
            traces = [RFTrace(trace=tr) for tr in stream.traces]
        #elif not traces is None and len(traces) > 0:
        #    traces = [RFTrace(trace=tr) for tr in traces]
        super(RFStream, self).__init__(traces=traces)

    def _write_test_header(self):
        for tr in self:
            tr._write_test_header()

    def write(self, filename, format, **kwargs):
        """
        Save stream to file including format specific headers.

        See :meth:`Stream.write() <obspy.core.stream.Stream.write>` in ObsPy.
        """
        for tr in self:
            tr._write_format_specific_header(format)
        super(RFStream, self).write(filename, format, **kwargs)

    def deconvolve(self, *args, **kwargs):
        """
        Deconvolve source component of stream from other components.

        All args and kwargs are passed to the function
        :func:`~rf.deconvolve.deconv`.
        """
        if len(self) % 3 != 0:
            raise ValueError('For deconvolution a 3 component stream is needed'
                             '. The provied stream is not divisible by 3.')
        i = 0
        while i < len(self):
            comps = ''.join(tr.stats.channel[-1] for tr in self[i:i + 3])
            if i == 0:
                comps0 = comps
            elif comps != comps0:
                raise ValueError('Error')
            deconv(self[i:i + 3], *args, **kwargs)
            i += 3

    def rf(self, method='P', filter=None, window=None, downsample=None,
           rotate='ZNE->LQT', deconvolve='time', **kwargs):
        r"""
        Calculate receiver functions in-place.

        :param method: 'P' for P receiver functions, 'S' for S receiver
            functions
        :param dictionary filter: filter stream with its
            :meth:`~obspy.core.stream.Stream.filter` method and given kwargs
        :type window: tuple of length 2
        :param window: trim stream relative to P- or S-onset
             with :meth:`~obspy.core.stream.Stream.trim` (seconds)
        :param float downsample: downsample stream with its
            :meth:`~obspy.core.stream.Stream.decimate` method
        :param rotate: 'ZNE->LQT' or 'ZNE->ZRT', rotate stream with its
            :meth:`~obspy.core.stream.Stream.rotate`
            method with the angles given by the back_azimuth and inclination
            attributes of the traces stats objects. You can set these to your
            needs or let them be computed by :func:`~rf.rfstream.rfstats`.
            The first component of the target component is assumed to
            be the source component for the deconvolution (L or Z).
        :param deconvolve: 'time' or 'freq' for time or frequency domain
            deconvolution by the streams
            :meth:`~rf.rfstream.RFStream.deconvolve`
            method. See :func:`~rf.deconvolve.deconv`,
            :func:`~rf.deconvolve.deconvt` and :func:`~rf.deconvolve.deconvf`
            for further documentation.
        :param \*\*kwargs: all other kwargs not mentioned here are
            passed to deconvolve

        After performing the deconvolution the Q/R and T components are
        multiplied by -1 to get a positive phase for a Moho-like positive
        velocity contrast. Furthermore for method='S' all components are
        mirrored at t=0 for a better comparison with P receiver functions.
        See source code of this function for the default
        deconvolution windows.
        """
        if method not in 'PS':
            raise NotImplementedError
        if filter:
            self.filter(**filter)
        if window:
            for tr in self:
                tr.trim(tr.stats.onset + window[0], tr.stats.onset + window[1])
        if downsample:
            for tr in self:
                if downsample <= tr.stats.sampling_rate:
                    tr.decimate(int(tr[0].stats.sampling_rate) // downsample)
        if isinstance(rotate, basestring):
            if rotate[-2:] == 'RT':
                self.rotate(rotate[1:5] + rotate[6:])
            else:
                self.rotate(rotate)
            src_comp = rotate.split('->')[-1][0]
        elif rotate:
            src_comp = rotate(self)
        if deconvolve:
            # set standard parameters for deconvolution
            if method == 'P' and deconvolve == 'time':
                kwargs.setdefault('winsrc', (-10, 30, 5))
                kwargs.setdefault('winrsp', (-20, 80))
                kwargs.setdefault('winrf', (-20, 80))
            elif method == 'S' and deconvolve == 'time':
                kwargs.setdefault('winsrc', (-10, 30, 5))
                kwargs.setdefault('winrsp', (-80, 20))
                kwargs.setdefault('winrf', (-80, 20))
            elif method == 'P':
                kwargs.setdefault('winsrc', (-20, 80, 5))
                kwargs.setdefault('tshift', 10)
            else:
                kwargs.setdefault('winsrc', (-20, 80, 5))
                kwargs.setdefault('tshift', 90)
            self.deconvolve(src_comp, method=deconvolve, **kwargs)
        for tr in self:
        # Mirrow Q/R and T component at 0s for S-receiver method for a better
        # comparison with P-receiver method (converted Sp wave arrives before
        # S wave, but converted Ps wave arrives after P wave)
            if method == 'S':
                tr.data = tr.data[::-1]
                tr.stats.onset = tr.stats.starttime + (tr.stats.endtime -
                                                       tr.stats.onset)
        # Multiply -1 on Q/R and T component, because Q/R component is pointing
        # towards the event after the rotation. For a positive phase at
        # a Moho-like velocity contrast, the Q/R component has to
        # point away from the event.
            if tr.stats.channel[-1] != src_comp:
                tr.data = -tr.data

    def moveout(self, phase='Ps', ref=6.4, model='iasp91'):
        """
        In-place moveout correction to a reference slowness.

        Needs stats attributes slowness and onset.

        :param phase: 'Ps', 'Sp', 'Ppss' or other multiples
        :param ref: reference ray parameter in s/deg
        :param model: path to model file or 'iasp91'
        """
        model = load_model(model)
        model.moveout(self, phase=phase, ref=ref)

    def ppoint(self, depth, phase='S', model='iasp91'):
        """
        Calculate coordinates of piercing point by 1D ray tracing.

        The iasp91 model is used. Piercing point coordinates are stored in the
        stats attributes plat and plon. Needs stats attributes
        station_latitude, station_longitude, slowness and back_azimuth.

        :param depth: depth of interface in km
        :param phase: 'P' for piercing points of P wave, 'S' for piercing
            points of S wave. Multiples are possible, too.
        :param model: path to model file or 'iasp91'
        :return: NumPy array with coordinates of piercing points

        .. note::

            `phase='S'` is usually wanted for P receiver functions and 'P'
            for S receiver functions.
        """
        model = load_model(model)
        for tr in self:
            model.ppoint(tr.stats, depth, phase=phase)
        return np.array([(tr.stats.plat, tr.stats.plon) for tr in self])

    def _moveout_xy(self, *args, **kwargs):
        for tr in self:
            tr._moveout_xy(*args, **kwargs)

    def _ppoint_xy(self, *args, **kwargs):
        for tr in self:
            tr._ppoint_xy(*args, **kwargs)


class RFTrace(Trace):
    """
    Class providing the Trace object for receiver function calculation.
    """
    def __init__(self, data=np.array([]), header={}, trace=None):
        if trace is not None:
            data = trace.data
            header = trace.stats
        super(RFTrace, self).__init__(data=data, header=header)
        self._read_format_specific_header()

    def _write_test_header(self):
        st = self.stats
        for head, val in zip(HEADERS, _HEADERS_EXAMPLE):
            if head in ('onset', 'event_time'):
                val = st.starttime + val
            st[head] = val

    def _read_format_specific_header(self, format=None):
        st = self.stats
        if format is None:
            if '_format' not in st:
                return
            format = st._format
        format = format.lower()
        if format == 'q':
            format = 'sh'
        try:
            header_map = zip(HEADERS, FORMATHEADERS[format])
        except KeyError:
            warnings.warn('Reading header of a file with this format is not '
                          'supported.')
            return
        for head, head_format in header_map:
            try:
                st[head] = st[format][head_format]
            except KeyError:
                continue
            try:
                convert = _HEADER_CONVERSIONS[format][head][0]
                st[head] = convert(st, head)
            except KeyError:
                pass

    def _write_format_specific_header(self, format):
        st = self.stats
        format = format.lower()
        if format == 'q':
            format = 'sh'
        try:
            header_map = zip(HEADERS, FORMATHEADERS[format])
        except KeyError:
            warnings.warn('Reading header of a file with this format is not '
                          'supported.')
            return
        if format not in st:
            st[format] = AttribDict({})
        for head, head_format in header_map:
            try:
                val = st[head]
            except KeyError:
                continue
            try:
                convert = _HEADER_CONVERSIONS[format][head][1]
                val = convert(st, head)
            except KeyError:
                pass
            st[format][head_format] = val

    def write(self, filename, format, **kwargs):
        """
        Save current trace into a file  including format specific headers.

        See :meth:`Trace.write() <obspy.core.trace.Trace.write>` in ObsPy.
        """
        RFStream([self]).write(filename, format, **kwargs)

    def _moveout_xy(self):
        """
        Depreciated! Moveout correction to a slowness of 6.4s/deg.

        The iasp91 model is used. The correction is independent from the type
        of receiver function. Needs stats attributes slowness and onset.
        """
        from rf import _xy
        st = self.stats
        print st
        self.data = _xy.psmout(self.data, st.slowness,
                               st.onset - st.starttime,
                               st.endtime - st.starttime, st.delta, 0)

    def _ppoint_xy(self, depth, method='P'):
        """
        Depreciated! Calculate coordinates of piercing point by 1D ray tracing.

        The iasp91 model is used. Piercing point coordinates are stored in the
        stats attributes `plat` and `plon`. Needs stats attributes
        station_latitude, station_longitude, slowness and back_azimuth.

        :param depth: depth of piercing points in km
        :param method: 'P' or 'S' for P or S waves
        """
        from rf import _xy
        if method not in 'PS':
            raise NotImplementedError()
        st = self.stats
        args = (depth, st.station_latitude, st.station_longitude,
                st.slowness, st.back_azimuth)
        pier_func = _xy.pspier if method == 'P' else _xy.sppier
        _, lat, lon = pier_func(*args)
        st.plat = lat
        st.plon = lon


def obj2stats(event=None, station=None):
    """
    Map event and station object to stats with attributes.

    :param event: ObsPy :class:`~obspy.core.event.Event` object
    :param station: station object with attributes latitude, longitude and
        elevation
    :return: ``stats`` object with station and event attributes
    """
    stats = AttribDict({})
    if event is not None:
        for key, getter in EVENT_GETTER:
            stats[key] = getter(event)
    if station is not None:
        for key, getter in STATION_GETTER:
            stats[key] = getter(station)
    return stats


def rfstats(stats=None, event=None, station=None, stream=None,
            phase='P', dist_range=None):
    """
    Calculate ray specific values like slowness for given event and station.

    :param stats: stats object with event and/or station attributes. Can be
        None if both event and station are given.
    :param event: ObsPy :class:`~obspy.core.event.Event` object
    :param station: station object with attributes latitude, longitude and
        elevation
    :param stream: If a stream is given, stats has to be None. In this case
        rfstats will be called for every stats object in the stream.
    :param phase: string with phase to look for in result of
        :func:`~obspy.taup.taup.getTravelTimes`. Usually this will be 'P' or
        'S' for P and S receiver functions, respectively.
    :type dist_range: tuple of length 2
    :param dist_range: if epicentral of event is not in this intervall, None
        is returned by this function,\n
        if phase == 'P' defaults to (30, 90),\n
        if phase == 'S' defaults to (50, 85)

    :return: ``stats`` object with event and station attributes, distance,
        back_azimuth, inclination, onset and slowness or None if epicentral
        distance is not in the given intervall
    """
    if stream is not None:
        assert stats is None
        for tr in stream:
            rfstats(tr.stats, event, station, None, phase, dist_range)
        return
    phase = phase.upper()
    if dist_range is None and phase in 'PS':
        dist_range = (30, 90) if phase == 'P' else (50, 85)
    if stats is None:
        stats = AttribDict({})
    stats.update(obj2stats(event=event, station=station))
    dist, baz, _ = gps2DistAzimuth(stats.station_latitude,
                                   stats.station_longitude,
                                   stats.event_latitude,
                                   stats.event_longitude)
    dist = kilometer2degrees(dist / 1000)
    if dist_range and not dist_range[0] <= dist <= dist_range[1]:
        return
    tts = getTravelTimes(dist, stats.event_depth)
    tts2 = getTravelTimes(dist, 0)
    tts = [tt for tt in tts if tt['phase_name'] == phase]
    tts2 = [tt for tt in tts2 if tt['phase_name'] == phase]
    if len(tts) == 0 or len(tts2) == 0:
        raise Exception('Taup does not return phase %s at event distance %s' %
                        (phase, dist))
    onset = stats.event_time + tts[0]['time']
    inc = tts2[0]['take-off angle']  # approximation
    v = 5.8 if 'P' in phase else 3.36  # iasp91
    slowness = 6371. * sin(pi / 180. * inc) / v / 180 * pi
    stats.update({'distance': dist, 'back_azimuth': baz, 'inclination': inc,
                  'onset': onset, 'slowness': slowness})
    return stats
