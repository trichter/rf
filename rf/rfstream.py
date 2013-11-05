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
try:
    from rf import _xy
except:
    warnings.warn("Error importing Fortran extensions module.")

STATION_GETTER = (('station_latitude', attrgetter('latitude')),
                  ('station_longitude', attrgetter('longitude')),
                  ('station_elevation', attrgetter('elevation')))

def __get(h1, h2): return lambda event: event[h1][0][h2]
EVENT_GETTER = (('event_latitude', __get('origins', 'latitude')),
                ('event_longitude', __get('origins', 'longitude')),
                ('event_depth', __get('origins', 'depth')),
                ('event_magnitude', __get('magnitudes', 'mag')),
                ('event_time', __get('origins', 'time')))
HEADERS = zip(*STATION_GETTER)[0] + zip(*EVENT_GETTER)[0] + (
              'onset', 'distance', 'back_azimuth', 'inclination', 'slowness')
FORMATHEADERS = {'sac': ('stla', 'stlo', 'stel', 'evla', 'evlo', 'evdp', 'mag',
                         'o', 'a', 'gcarc', 'baz', 'user0', 'user1'),
                 # fields 'dcvreg', 'dcvinci' and 'pwdw' are violated for
                 # station information
                 'sh': ('DCVREG', 'DCVINCI', 'PWDW', 'LAT', 'LON', 'DEPTH',
                        'MAGNITUDE', 'ORIGIN', 'P-ONSET', 'DISTANCE',
                        'AZIMUTH', 'INCI', 'SLOWNESS')}

def __rel2UTC(stats, head):
    return stats.starttime + stats[head]
def __UTC2rel(stats, head):
    return stats[head] - stats.starttime
_HEADER_CONVERSIONS = {'sac': {'onset': (__rel2UTC, __UTC2rel),
                               'event_time': (__rel2UTC, __UTC2rel)}}


class RFStream(Stream):
    """
    Class providing the necessary functions for receiver function calculation.

    To initialize a RFStream from a ObsPy stream use

    >>> rfstream = RFStream(stream=obspy_stream)
    """
    def __init__(self, traces=None, stream=None):
        if stream is not None:
            traces = stream.traces
        if not traces is None and len(traces) > 0:
            traces = [RFTrace(trace=tr) for tr in traces]
        super(RFStream, self).__init__(traces=traces)

    def write(self, filename, format, **kwargs):
        """
        Saves stream to file including format specific headers.

        See :meth:`Stream.write() <obspy.core.stream.Stream.write>` in ObsPy.
        """
        # Check all traces for masked arrays and raise exception.
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
            comps = ''.join(tr.stats.component for tr in self[i, i + 3])
            if i == 0:
                comps0 = comps
            elif comps != comps0:
                raise ValueError('Error')
            deconv(self[i, i + 3], *args, **kwargs)
            i += 3

    def rf(self, method='P', window=None, downsample=None,
           filter=None,
           rotate='ZNE->LQT', deconvolve='time', **deconvolve_kwargs):
        """
        TODO: DOC
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
        if rotate:
            self.rotate(rotate)
            src_comp = rotate.split('->')[-1][0]
        if deconvolve:
            # set standard parameters for deconvolution
            kwargs = deconvolve_kwargs
            if method == 'P' and deconvolve == 'time':
                kwargs.set_default('winsrc', (-10, 30, 5))
                kwargs.set_default('winrsp', (-20, 80))
                kwargs.set_default('winrf', (-20, 80))
            elif method == 'S' and deconvolve == 'time':
                kwargs.set_default('winsrc', (-10, 30, 5))
                kwargs.set_default('winrsp', (-80, 20))
                kwargs.set_default('winrf', (-80, 20))
            elif method == 'P':
                kwargs.set_default('winsrc', (-20, 80, 5))
                kwargs.set_default('tshift', 10)
            else:
                kwargs.set_default('winsrc', (-20, 80, 5))
                kwargs.set_default('tshift', 90)
            self.deconvolve(src_comp, method=deconvolve, **deconvolve_kwargs)
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


    def moveout(self):
        """
        Moveout correction to a slowness of 6.4s/deg.

        The iasp91 model is used.
        """
        for tr in self:
            tr.moveout()

    def ppoint(self, depth, method='P'):
        """
        Calculate coordinates of piercing point at `depth` by 1D ray tracing.

        The iasp91 model is used. Piercing point coordinates are stored in the
        stats attributes `plat` and `plon`.
        """
        for tr in self:
            tr.piercing_points(depth, method=method)


class RFTrace(Trace):
    def __init__(self, data=np.array([]), header={}, trace=None):
        if trace is not None:
            data = trace.data
            header = trace.stats
        super(RFTrace, self).__init__(data=data, header=header)
        self._read_format_specific_header()

    def write(self, filename, format, **kwargs):
        RFStream([self]).write(filename, format, **kwargs)

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

    def moveout(self):
        st = self.stats
        self.data = _xy.psmout(self.data, st.slowness,
                                   st.ponset - st.starttime,
                                   st.endtime - st.starttime, st.dt, 0)

    def ppoint(self, depth, method='P'):
        if method not in 'PS':
            raise NotImplementedError()
        st = self.stats
        args = (depth, st.station.latitude, st.station.longitude,
                st.slowness, st.baz)
        pier_func = _xy.pspier if method == 'P' else _xy.sppier
        _, lat, lon = pier_func(*args)
        st.plat = lat
        st.plon = lon


def obj2stats(event=None, station=None):
    stats = AttribDict({})
    if event is not None:
        for key, getter in EVENT_GETTER:
            stats[key] = getter(event)
    if station is not None:
        for key, getter in STATION_GETTER:
            stats[key] = getter(station)
    return stats

def rfstats(stats=None, event=None, station=None, phase='P', dist_range=None):
    """
    Calculate important rf-specific values.

    TODO: DOC
    """
    phase = phase.upper()
    if dist_range is None:
        dist_range = (30, 90) if 'P' in phase else (50, 85)
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
