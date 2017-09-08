# Copyright 2013-2016 Tom Eulenfeld, MIT license
"""
Utility functions and classes for receiver function calculation.
"""
import collections
import inspect
import itertools
from pkg_resources import resource_filename

from decorator import decorator
import numpy as np


DEG2KM = 111.2  #: Conversion factor from degrees epicentral distance to km


def _get_stations(inventory):
    channels = inventory.get_contents()['channels']
    stations = {ch[:-1] + '?': ch[-1] for ch in channels}
    return stations


def iter_event_data(events, inventory, get_waveforms, phase='P',
                    request_window=None, pad=10, pbar=None, **kwargs):
    """
    Return iterator yielding three component streams per station and event.

    :param events: list of events or `~obspy.core.event.Catalog` instance
    :param inventory: `~obspy.core.inventory.inventory.Inventory` instance
        with station and channel information
    :param get_waveforms: Function returning the data. It has to take the
        arguments network, station, location, channel, starttime, endtime.
    :param phase: Considered phase, e.g. 'P', 'S', 'PP'
    :type request_window: tuple (start, end)
    :param request_window: requested time window around the onset of the phase
    :param float pad: add specified time in seconds to request window and
       trim afterwards again
    :param pbar: tqdm_ instance for displaying a progressbar
    :param kwargs: all other kwargs are passed to `~rf.rfstream.rfstats()`

    :return: three component streams with raw data

    Example usage with progressbar::

        from tqdm import tqdm
        from rf.util import iter_event_data
        with tqdm() as t:
            for stream3c in iter_event_data(*args, pbar=t):
                do_something(stream3c)

    .. _tqdm: https://pypi.python.org/pypi/tqdm
    """
    from rf.rfstream import rfstats, RFStream
    method = phase[-1].upper()
    if request_window is None:
        request_window = (-50, 150) if method == 'P' else (-100, 50)
    stations = _get_stations(inventory)
    if pbar is not None:
        pbar.total = len(events) * len(stations)
    for event, seedid in itertools.product(events, stations):
        if pbar is not None:
            pbar.update(1)
        origin_time = (event.preferred_origin() or event.origins[0])['time']
        try:
            args = (seedid[:-1] + stations[seedid], origin_time)
            coords = inventory.get_coordinates(*args)
        except:  # station not available at that time
            continue
        stats = rfstats(station=coords, event=event, phase=phase, **kwargs)
        if not stats:
            continue
        net, sta, loc, cha = seedid.split('.')
        starttime = stats.onset + request_window[0]
        endtime = stats.onset + request_window[1]
        kws = {'network': net, 'station': sta, 'location': loc,
               'channel': cha, 'starttime': starttime - pad,
               'endtime': endtime + pad}
        try:
            stream = get_waveforms(**kws)
        except:  # no data available
            continue
        stream.trim(starttime, endtime)
        stream.merge()
        if len(stream) != 3:
            from warnings import warn
            warn('Need 3 component seismograms. %d components '
                 'detected for event %s, station %s.'
                 % (len(stream), event.resource_id, seedid))
            continue
        if any(isinstance(tr.data, np.ma.masked_array) for tr in stream):
            from warnings import warn
            warn('Gaps or overlaps detected for event %s, station %s.'
                 % (event.resource_id, seedid))
            continue
        for tr in stream:
            tr.stats.update(stats)
        yield RFStream(stream)


def iter_event_metadata(events, inventory, pbar=None):
    """
    Return iterator yielding metadata per station and event.

    :param events: list of events or `~obspy.core.event.Catalog` instance
    :param inventory: `~obspy.core.inventory.inventory.Inventory` instance
        with station and channel information
    :param pbar: tqdm_ instance for displaying a progressbar
    """
    stations = _get_stations(inventory)
    if events is None:
        events = [None]
    if pbar is not None:
        pbar.total = len(events) * len(stations)
    for event, seedid in itertools.product(events, stations):
        if pbar is not None:
            pbar.update(1)
        net, sta, loc, cha = seedid.split('.')
        meta = {'network': net, 'station': sta, 'location': loc,
                'channel': cha}
        if event is not None:
            ot = (event.preferred_origin() or event.origins[0])['time']
            meta['event_time'] = ot
        yield meta


class IterMultipleComponents(object):

    """
    Return iterable to iterate over associated components of a stream.

    :param stream: Stream with different, possibly many traces. It is
        split into substreams with the same seed id (only last character
        i.e. component may vary)
    :type key: str or None
    :param key: Additionally, the stream is grouped by the values of
         the given stats entry to differentiate between e.g. different events
         (for example key='starttime', key='onset')
    :type number_components: int, tuple of ints or None
    :param number_components: Only iterate through substreams with
         matching number of components.
    """

    def __init__(self, stream, key=None, number_components=None):
        substreams = collections.defaultdict(stream.__class__)
        for tr in stream:
            k = (tr.id[:-1], str(tr.stats[key]) if key is not None else None)
            substreams[k].append(tr)
        n = number_components
        self.substreams = [s for _, s in sorted(substreams.items())
                           if n is None or len(s) == n or len(s) in n]

    def __len__(self):
        return len(self.substreams)

    def __iter__(self):
        for s in self.substreams:
            yield s


def direct_geodetic(latlon, azi, dist):
    """
    Solve direct geodetic problem with geographiclib.

    :param tuple latlon: coordinates of first point
    :param azi: azimuth of direction
    :param dist: distance in km

    :return: coordinates (lat, lon) of second point on a WGS84 globe
    """
    from geographiclib.geodesic import Geodesic
    coords = Geodesic.WGS84.Direct(latlon[0], latlon[1], azi, dist * 1000)
    return coords['lat2'], coords['lon2']


__CACHE = {}


def minimal_example_rf():
    """
    Return receiver functions calculated from the data returned by read_rf().
    """
    cache_key = 'minimal_example_rf'
    if cache_key in __CACHE:
        return __CACHE[cache_key].copy()
    from rf.rfstream import read_rf, rfstats
    stream = read_rf()
    rfstats(stream)
    stream.filter('bandpass', freqmin=0.5, freqmax=2)
    stream.trim2(10, 110, reftime='starttime')
    stream.rf(winsrc=(-5, 25, 5))
    stream.moveout()
    stream.trim2(-10, 80, reftime='onset')
    stream.ppoints(50)
    __CACHE[cache_key] = stream
    return stream.copy()


def minimal_example_Srf():
    """
    Return S receiver functions calculated from example data.
    """
    cache_key = 'minimal_example_Srf'
    if cache_key in __CACHE:
        return __CACHE[cache_key].copy()
    from rf.rfstream import read_rf, rfstats
    fname = resource_filename('rf', 'example/minimal_example_S.tar.gz')
    stream = read_rf(fname)
    rfstats(stream, phase='S')
    stream.filter('bandpass', freqmin=0.2, freqmax=0.5)
    stream.trim2(10, 120, reftime='starttime')
    stream.rf(method='S', winsrc=(-5, 15, 5))
    stream.moveout(phase='Sp')
    stream.ppoints(50, pp_phase='P')
    __CACHE[cache_key] = stream
    return stream.copy()


@decorator
def _add_processing_info(func, *args, **kwargs):
    from rf import __version__
    args_ = inspect.getcallargs(func, *args, **kwargs)
    if args_.pop('self', None) is None:
        args_.pop('stream')
    kw = args_.pop('kwargs', {})
    kw.update(args_)
    boxes = kw.pop('boxes', None)
    if boxes:
        kw['latlon'] = boxes[0]['profile']['latlon']
        kw['azimuth'] = boxes[0]['profile']['azimuth']
        kw['length'] = boxes[0]['profile']['length']
    info = 'rf {version}: {function}(%s)'.format(
        version=__version__, function=func.__name__)
    arguments = ['%s=%s' % (k, repr(v)) if not isinstance(v, str) else
                 "%s='%s'" % (k, v) for k, v in kw.items()]
    info = info % '::'.join(sorted(arguments))
    stream = func(*args, **kwargs)
    try:
        for tr in stream:
            tr._internal_add_processing_info(info)
    except:
        pass
    return stream
