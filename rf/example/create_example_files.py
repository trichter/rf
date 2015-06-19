from obspy import readEvents, read_inventory
from obspy.core import Stream, UTCDateTime as UTC
from obspy.fdsn import Client
from obspy.arclink import Client as ArcClient

evname = './example_events.xml'
invname = './example_inventory.xml'
wavname = './example_data.mseed'
wavformat = 'mseed'

lon, lat = -70, -21
t1, t2 = UTC('2011-01-01'), UTC('2011-06-01')
seed_id = 'CX.PB01..BH?'


def get_events():
    try:
        return readEvents(evname)
    except:
        pass
    client = Client()
    events = client.get_events(starttime=t1, endtime=t2, latitude=lat,
                               longitude=lon, minradius=20, maxradius=100,
                               minmagnitude=6.)
    events.write(evname, 'QUAKEML')
    return events


def get_inventory():
    try:
        return read_inventory(invname)
    except:
        pass
    client = Client('GFZ')
    net, sta, loc, cha = seed_id.split('.')
    inv = client.get_stations(starttime=t1, endtime=t2, network=net,
                              station=sta, location=loc, channel=cha,
                              level='channel')
                              # latitude=lat, longitude=lon, maxradius=10)
    inv.write(invname, 'STATIONXML')
    return inv


def get_waveforms():
    events = get_events()
    client = ArcClient()
    wforms = Stream()
    for event in events:
        t = event.preferred_origin().time
        args = seed_id.split('.') + [t + 5 * 60, t + 14 * 60]
        wforms.extend(client.getWaveform(*args))
    wforms.decimate(int(round(wforms[0].stats.sampling_rate)) // 5,
                    no_filter=True)
    wforms.write(wavname, wavformat)

print get_inventory()
print get_events()
get_waveforms()