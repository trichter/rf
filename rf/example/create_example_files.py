# Copyright 2013-2016 Tom Eulenfeld, MIT license
from obspy import read_events, read_inventory, Stream, UTCDateTime as UTC
from obspy.clients.fdsn import Client

from rf.rfstream import obj2stats, rfstats
from rf import RFStream


evname = './example_events.xml'
invname = './example_inventory.xml'
wavname = './example_data.mseed'
wavname2 = './minimal_example.sac'

lon, lat = -70, -21
t1, t2 = UTC('2011-02-01'), UTC('2011-06-01')
seedid = 'CX.PB01..BH?'


def get_events():
    try:
        return read_events(evname)
    except Exception:
        pass
    client = Client()
    events = client.get_events(starttime=t1, endtime=t2, latitude=lat,
                               longitude=lon, minradius=30, maxradius=90,
                               minmagnitude=6., maxmagnitude=6.5)
    events.write(evname, 'QUAKEML')
    return events


def get_inventory():
    try:
        return read_inventory(invname)
    except Exception:
        pass
    client = Client('GFZ')
    net, sta, loc, cha = seedid.split('.')
    inv = client.get_stations(starttime=t1, endtime=t2, network=net,
                              station=sta, location=loc, channel=cha,
                              level='channel')
#                               latitude=lat, longitude=lon, maxradius=10)
    inv.write(invname, 'STATIONXML')
    return inv


def get_waveforms():
    events = get_events()[::-1]
    client = Client('GFZ')
    stream_raw = Stream()
    stream = RFStream()
    coords = inventory.get_coordinates(seedid[:-1] + 'Z')
    for i, event in enumerate(events):
        t = event.preferred_origin().time
        args = seedid.split('.') + [t + 4.9 * 60, t + 14.1 * 60]
        s = client.get_waveforms(*args)
        s.trim(t+5*60, t+14*60)
        s.decimate(int(round(s[0].stats.sampling_rate)) // 5, no_filter=True)
        stream_raw.extend(s)
        if i in (0, 2, 4):
            s = s.copy()
            stats = rfstats(station=coords, event=event, dist_range=(20, 95))
            if stats is None:
                continue
            s.trim(stats.onset - 25, stats.onset + 75)
            stats = obj2stats(station=coords, event=event)
            s = RFStream(s)
            for tr in s:
                tr.stats.update(stats)
            stream.extend(s)
    stream_raw.write(wavname, 'MSEED')
    stream.write(wavname2, 'SAC')

inventory = get_inventory()
events = get_events()
get_waveforms()
