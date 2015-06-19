# -*- coding: utf-8 -*-
"""
Classes and functions for receiver function calculation.
"""

import logging
from operator import itemgetter
import warnings

import numpy as np
from obspy import read
from obspy.core import Stream, Trace
from obspy.core.util import AttribDict
from obspy.core.util.geodetics import gps2DistAzimuth, kilometer2degrees
from obspy.taup import TauPyModel
from rf.deconvolve import deconv
from rf.simple_model import load_model

def __get_event_origin(h):
    return lambda event: event.preferred_origin()[h]


def __rel2UTC(stats, head):
    return stats.starttime + stats[head]


def __UTC2rel(stats, head):
    return stats[head] - stats.starttime


STATION_GETTER = (('station_latitude', itemgetter('latitude')),
                  ('station_longitude', itemgetter('longitude')),
                  ('station_elevation', itemgetter('elevation')))
EVENT_GETTER = (
    #('event_id', lambda event: _get_event_id(event)),
    ('event_latitude', __get_event_origin('latitude')),
    ('event_longitude', __get_event_origin('longitude')),
    ('event_depth', lambda event: event.preferred_origin()['depth'] / 1000.),
    ('event_magnitude', lambda event: event.preferred_magnitude()['mag']),
    ('event_time', __get_event_origin('time')))
HEADERS = zip(*STATION_GETTER)[0] + zip(*EVENT_GETTER)[0] + (
    'onset', 'distance', 'back_azimuth', 'inclination', 'slowness')
FORMATHEADERS = {'sac': ('stla', 'stlo', 'stel', 'evla', 'evlo',
                         'evdp', 'mag',
                         'o', 'a', 'gcarc', 'baz', 'user0', 'user1'),
                 # fields 'dcvreg', 'dcvinci' and 'pwdw' are violated for
                 # station information
                 'sh': ('DCVREG', 'DCVINCI', 'PWDW',
                        'LAT', 'LON', 'DEPTH',
                        'MAGNITUDE', 'ORIGIN', 'P-ONSET', 'DISTANCE',
                        'AZIMUTH', 'INCI', 'SLOWNESS')}
_HEADER_CONVERSIONS = {'sac': {'onset': (__rel2UTC, __UTC2rel),
                               'event_time': (__rel2UTC, __UTC2rel)}}
_HEADERS_EXAMPLE = (50.3, -100.2, 400.3,
                    -20.32, 10., 12.4, 6.5, -40.432,
                    20.643, 57.6, 90.1, 10.2, 10.)

_TF = '.datetime:%Y-%m-%dT%H:%M:%S'
H5INDEX = ('{network}.{station}.{location}/{event_time%s}/' % _TF +
           '{channel}_{starttime%s}_{endtime%s}' % (_TF, _TF))
H5INDEX_STACK = '{network}.{station}.{location}/{channel}'


def set_index(index='rf'):
    import obspyh5
    if index == 'rf':
        index = H5INDEX
    elif index == 'rf_stack':
        index = H5INDEX_STACK
    obspyh5.set_index(index)


def read_rf(fname=None, format_=None, **kwargs):
    """
    Read waveform files into RFStream object.

    See :func:`read() <obspy.core.stream.read>` in ObsPy.
    """
    return RFStream(read(fname, format=format_, **kwargs))


class RFStream(Stream):

    """
    Class providing the necessary functions for receiver function calculation.

    To initialize a RFStream from a Stream object use

    >>> rfstream = RFStream(stream)

    To initialize a RFStream from a file use

    >>> rfstream = read_rf('test.SAC')

    Format specific headers are loaded into the stats object of all traces.
    """

    def __init__(self, traces=None, warn=True):
        self.traces = []
        if isinstance(traces, Trace):
            traces = [traces]
        if traces:
            for tr in traces:
                if not isinstance(tr, RFTrace):
                    tr = RFTrace(trace=tr, warn=warn)
                self.traces.append(tr)

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
            if format.upper() == 'Q':
                tr.stats.station = tr.id
        super(RFStream, self).write(filename, format, **kwargs)
        if format.upper() == 'Q':
            for tr in self:
                tr.stats.station = tr.stats.station.split('.')[1]

    def deconvolve(self, method='P', deconvolve_method='time', **kwargs):
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
            # set standard parameters for deconvolution
            stats = self[i].stats
            lensec = stats.endtime - stats.starttime
            onset = stats.onset - stats.starttime
            if method == 'P' and deconvolve_method == 'time':
                def_kwargs = {'winsrc': (-10, 30, 5),
                              'winrsp': (-onset, lensec - onset),
                              'winrf': (-onset, lensec - onset)}
            elif method == 'S' and deconvolve_method == 'time':
                def_kwargs = {'winsrc': (-10, 30, 5),
                              'winrsp': (onset - lensec, onset),
                              'winrf': (onset - lensec, onset)}
            elif method == 'P':
                def_kwargs = {'winsrc': (-onset + 5, lensec - onset - 5, 5),
                              'tshift': onset}
            else:
                def_kwargs = {'winsrc': (onset - lensec + 5, onset - 5),
                              'tshift': lensec - onset}
            nkwargs = kwargs.copy()
            for k in def_kwargs:
                if k not in kwargs:
                    nkwargs[k] = def_kwargs[k]
            try:
                deconv(self[i:i + 3], method=deconvolve_method, **nkwargs)
            except:
                print('error while calculating the deconvolution')
                for tr in self[i:i + 3]:
                    self.remove(tr)
                continue
            i += 3

    def rf(self, method='P', filter=None, window=None, downsample=None,
           rotate='ZNE->LQT', source_component='L',
           deconvolve='time', **kwargs):
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
        :param rotate: 'ZNE->LQT' or 'NE->RT', rotate stream with its
            :meth:`~obspy.core.stream.Stream.rotate`
            method with the angles given by the back_azimuth and inclination
            attributes of the traces stats objects. You can set these to your
            needs or let them be computed by :func:`~rf.rfstream.rfstats`.
        :param source_component: character of the source component
            (i.e. 'L' or 'Z' depending on rotaion)
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
                    tr.decimate(int(tr.stats.sampling_rate) // downsample)
        if isinstance(rotate, basestring):
            self.rotate(rotate)
        elif rotate:
            rotate(self)
        if deconvolve:
            self.deconvolve(method=method, deconvolve_method=deconvolve,
                            source_component=source_component, **kwargs)
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
            if tr.stats.channel[-1] != source_component:
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

    def stack(self):
        """
        Stack traces with the same id.

        Traces with same id need to have the same number of datapoints.
        """
        ids = set(tr.id for tr in self)
        traces = []
        for id in ids:
            net, sta, loc, cha = id.split('.')
            data = np.mean([tr.data for tr in self if tr.id == id], axis=0)
            header = {'network': net, 'station': sta, 'location': loc,
                      'channel': cha, 'sampling_rate': tr.stats.sampling_rate}
            onset = tr.stats.onset - tr.stats.starttime
            tr = RFTrace(data=data, header=header)
            tr.stats['onset'] = tr.stats['starttime'] + onset
            traces.append(tr)
        self.traces = traces

    def plot_rf(self, fname=None, norm=1., fig_width=7., trace_height=0.5,
                stack_height=0.5,
                fill=False, window=None, downsample=None, title=True,
                info=[('back_azimuth', u'baz (°)', 'b'),
                      ('distance', u'dist (°)', 'r')]):
        """
        Create receiver function plot.

        :param fname: Filename to save plot to. Can be None. In this case
            the figure is left open.
        :param fig_width: Width of figure in inches.
        :param trace_height: Height of one trace in inches.
        :param fill: Waether to fill receiver functions or not.
        :param downsample: Downsample to frequency (in Hz) with
            Stream.decimate. Filtering is not performed. When saving in a
            vector format the plot size can be reduced in this way.
        :param title: Print seed id as a title.
        :param info: Plot one additional axes showing maximal two entries of
            the stats object. Each entry in this list is a list consisting of
            three entries: key, label and color.
            info can be None. In this case no additional axes is plotted.
        """
        import matplotlib.pyplot as plt
        from matplotlib.ticker import AutoMinorLocator, MaxNLocator
        if len(self) == 0:
            return
        if window:
            for tr in self:
                tr.trim(tr.stats.onset + window[0], tr.stats.onset + window[1])
        if downsample:
            for tr in self:
                tr.decimate(int(round(tr.stats.sampling_rate)) // downsample,
                            no_filter=True)
        # calculate lag times
        stats = self[0].stats
        N = len(self)
        t0 = stats.onset - stats.starttime
        t2 = stats.endtime - stats.starttime
        times = np.linspace(-t0, t2 - t0, stats.npts, endpoint=True)
        # calculate axes and figure dimensions
        # big letters: inches, small letters: figure fraction
        H = trace_height
        HS = stack_height
        FB = 0.5
        FT = 0.2
        DW = 0.1
        FH = H * (N + 2) + HS + FB + FT + DW
        h = H / FH
        hs = HS / FH
        fb = FB / FH
        ft = FT / FH
        FL = 0.5
        FR = 0.2
        FW = fig_width
        FW3 = 0.8
        FW2 = FW - FL - FR - (DW + FW3) * bool(info)
        fl = FL / FW
        fr = FR / FW
        fw2 = FW2 / FW
        fw3 = FW3 / FW
        # init figure and axes
        fig = plt.figure(figsize=(FW, FH))
        ax1 = fig.add_axes([fl, fb, fw2, h * (N + 2)])
        ax2 = fig.add_axes([fl, 1 - ft - hs, fw2, hs], sharex=ax1)
        if info:
            ax3 = fig.add_axes(
                [1 - fr - fw3, fb, fw3, h * (N + 2)], sharey=ax1)
            info = list(info)
            info[0] = [ax3] + list(info[0])
            if len(info) > 1:
                ax4 = ax3.twiny()
                info[1] = [ax4] + list(info[1])
        # plot stack and individual receiver functions
        stack = self.copy()
        stack.stack()
        if len(stack) > 1:
            warnings.warn('Different stations in one RF plot.')

        def _rf_fill(ax, t, d, i):
            ax.fill_between(t, d + i, i, where=d >= 0, lw=0., facecolor='k')
            ax.fill_between(t, d + i, i, where=d < 0, lw=0, facecolor='grey')

        def _plot(ax, t, d, i):
            if fill:
                _rf_fill(ax, t, d, i)
                ax.plot(t, d + i, 'k')
            else:
                ax.plot(t, d + i, 'k')
        _plot(ax2, times, stack[0].data, 0)
        for i, tr in enumerate(self):
            _plot(ax1, times, tr.data * norm, i + 1)
        # plot right axes with header information
        for ax, header, label, color in info:
            data = [tr.stats[header] for tr in self]
            ax.plot(data, 1 + np.arange(len(self)), '.' + color, mec=color)
            ax.set_xlabel(label, color=color, size='small')
            if header == 'back_azimuth':
                ax.set_xticks(np.arange(5) * 90)
                ax.set_xticklabels(['0', '', '180', '', '360'], size='small')
            else:
                ax.xaxis.set_major_locator(MaxNLocator(4))
                for l in ax.get_xticklabels():
                    l.set_fontsize('small')
            ax.xaxis.set_minor_locator(AutoMinorLocator())
        # set x and y limits
        ax1.set_xlim(times[0], times[-1])
        ax1.set_ylim(-0.5, N + 1.5)
        ax1.set_yticklabels('')
        ax1.set_xlabel('time (s)')
        ax1.xaxis.set_minor_locator(AutoMinorLocator())
        for l in ax2.get_xticklabels():
            l.set_visible(False)
        ax2.yaxis.set_major_locator(MaxNLocator(4))
        for l in ax2.get_yticklabels():
            l.set_fontsize('small')
        # plot title and save plot
        if title:
            bbox = dict(boxstyle='round', facecolor='white', alpha=0.8, lw=0)
            text = '%s traces  %s' % (len(self), stack[0].id)
            ax2.annotate(text, (1 - 0.5 * fr, 1 - 0.5 * ft),
                         xycoords='figure fraction', va='top', ha='right',
                         bbox=bbox, clip_on=False)
        if fname:
            fig.savefig(fname)
            plt.close(fig)

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

    def __init__(self, data=np.array([]), header={}, trace=None, warn=True):
        if trace is not None:
            data = trace.data
            header = trace.stats
        super(RFTrace, self).__init__(data=data, header=header)
        st = self.stats
        if ('_format'in st and st._format.upper() == 'Q' and
            st.station.count('.') > 0):
                st.network, st.station, st.location = st.station.split('.')[:3]
        self._read_format_specific_header(warn=warn)

    def __str__(self, id_length=None):
        out = (u' | {event_magnitude:.1f}M dist:{distance:.1f} '
               u'baz:{back_azimuth:.1f}')
        try:
            out = out.format(**self.stats)
        except KeyError:
            out = ''
        return super(RFTrace, self).__str__(id_length=id_length) + out

    def _write_test_header(self):
        st = self.stats
        for head, val in zip(HEADERS, _HEADERS_EXAMPLE):
            if head in ('onset', 'event_time'):
                val = st.starttime + val
            st[head] = val

    def _read_format_specific_header(self, format=None, warn=True):
        st = self.stats
        if format is None:
            if '_format' not in st:
                return
            format = st._format
        format = format.lower()
        if format == 'q':
            format = 'sh'
        if format == 'h5':
            return
        try:
            header_map = zip(HEADERS, FORMATHEADERS[format])
        except KeyError:
            if warn:
                warnings.warn('Reading rf header of a file with this format '
                              'is not supported.')
            return
        for head, head_format in header_map:
            try:
                value = st[format][head_format]
            except KeyError:
                continue
            else:
                if format == 'sac' and '-12345' in str(value):
                    pass
                else:
                    st[head] = value
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
        if format == 'h5':
            return
        try:
            header_map = zip(HEADERS, FORMATHEADERS[format])
        except KeyError:
            if format != 'h5':
                msg = ("rf in-/output of file format '%s' is not supported" %
                       format)
                warnings.warn(msg)
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
            phase='P', dist_range=None, model='iasp91'):
    """
    Calculate ray specific values like slowness for given event and station.

    :param stats: stats object with event and/or station attributes. Can be
        None if both event and station are given.
    :param event: ObsPy :class:`~obspy.core.event.Event` object
    :param station: station object with attributes latitude, longitude and
        elevation
    :param stream: If a stream is given, stats has to be None. In this case
        rfstats will be called for every stats object in the stream.
    :param phase: string with phase. Usually this will be 'P' or
        'S' for P and S receiver functions, respectively.
    :type dist_range: tuple of length 2
    :param dist_range: if epicentral of event is not in this intervall, None
        is returned by this function,\n
        if phase == 'P' defaults to (30, 90),\n
        if phase == 'S' defaults to (50, 85)
    :param model: model for travel time calculation. See the :mod:`obspy.taup`
        module.
    :return: ``stats`` object with event and station attributes, distance,
        back_azimuth, inclination, onset and slowness or None if epicentral
        distance is not in the given intervall
    """
    if stream is not None:
        assert stats is None
        for tr in stream:
            rfstats(tr.stats, event, station, None, phase, dist_range, model)
        return
    phase = phase.upper()
    if dist_range is None and phase in 'PS':
        dist_range = (30, 90) if phase == 'P' else (50, 85)
    if stats is None:
        stats = AttribDict({})
    if event is not None and station is not None:
        stats.update(obj2stats(event=event, station=station))
    dist, baz, _ = gps2DistAzimuth(stats.station_latitude,
                                   stats.station_longitude,
                                   stats.event_latitude,
                                   stats.event_longitude)
    dist = kilometer2degrees(dist / 1000)
    if dist_range and not dist_range[0] <= dist <= dist_range[1]:
        return
    model = TauPyModel(model=model)
    arrivals = model.get_travel_times(stats.event_depth, dist, (phase,))
    if len(arrivals) == 0:
        raise Exception('TauPy does not return phase %s at distance %s' %
                        (phase, dist))
    if len(arrivals) > 1:
        from warnings import warn
        msg = ('TauPy returns more than one arrival for phase %s at '
               'distance -> take first arrival' )
        warn(msg % (phase, dist))
    arrival = arrivals[0]
    onset = stats.event_time + arrival.time
    inc = arrival.incident_angle
    slowness = arrival.ray_param_sec_degree
    stats.update({'distance': dist, 'back_azimuth': baz, 'inclination': inc,
                  'onset': onset, 'slowness': slowness})
    return stats
