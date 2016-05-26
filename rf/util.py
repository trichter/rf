import collections
import itertools


DEG2KM = 111.2


class IterEventData(object):

    def __init__(self, catalog, inventory, get_waveforms, phase='P',
                 request_window=None, pad=10, **kwargs):
        self.catalog = catalog
        self.inventory = inventory
        self.get_waveforms = get_waveforms
        method = phase[-1].upper()
        if request_window is None:
            request_window = (-50, 150) if method == 'P' else (-100, 50)
        self.request_window = request_window
        self.phase = phase
        self.pad = pad
        self.kwargs = kwargs
        channels = inventory.get_contents()['channels']
        self.stations = {ch[:-1] + '?': ch[-1] for ch in channels}

    def __len__(self):
        return len(self.catalog) * len(self.stations)

    def __iter__(self):
        from rf import rfstats, RFStream
        for event, seedid in itertools.product(self.catalog, self.stations):
            origin_time = (event.preferred_origin() or
                           event.origins[0])['time']
            try:
                gc = self.inventory.get_coordinates
                coords = gc(seedid[:-1] + self.stations[seedid], origin_time)
            except:  # station not availlable at that time
                # todo log
                continue
            stats = rfstats(station=coords, event=event,
                            phase=self.phase, **self.kwargs)
            if not stats:
                continue
            net, sta, loc, cha = seedid.split('.')
            starttime = stats.onset + self.request_window[0]
            endtime = stats.onset + self.request_window[1]
            kws = {'network': net, 'station': sta, 'location': loc,
                   'channel': cha, 'starttime': starttime - self.pad,
                   'endtime': endtime + self.pad}
            try:
                stream = self.get_waveforms(**kws)
            except:  # no data availlable
                # todo log
                continue
            stream.trim(starttime, endtime)
            stream.merge()
            if len(stream) != 3:
                # todo log
                from warnings import warn
                warn('Need 3 component seismograms. %d components '
                     'detected for event %s, station %s.'
                     % (len(stream), event.resource_id, seedid))
                continue
            for tr in stream:
                tr.stats.update(stats)
            yield RFStream(stream, warn=False)


class IterMultipleComponents(object):

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


def direct_geodetic(lonlat, azi, dist):
    from geographiclib.geodesic import Geodesic
    coords = Geodesic.WGS84.Direct(lonlat[1], lonlat[0], azi, dist * 1000)
    return coords['lon2'], coords['lat2']
