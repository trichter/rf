# by TR

from obspy import read
from obspy.core.util import AttribDict
from obspy.core.util.geodetics import gps2DistAzimuth, kilometer2degrees
from obspy.taup.taup import getTravelTimes
from rf import io, conf
from rf.deconvolve import deconv
import glob
import os

io.set_paths('~/obspyDMT-data/2010-01-01_2010-01-10_5.5_9.9')


def event2stats(lat, lon, event, phase='P', dist_range=(30, 90)):
    phase = phase.upper()
    ori = event.origins[0]
    dist, baz, _ = gps2DistAzimuth(lat, lon,
                                   ori.latitude, ori.longitude)
    dist = kilometer2degrees(dist / 1000)
    if not dist_range[0] <= dist <= dist_range[1]:
        return
    tts = getTravelTimes(dist, ori.depth)
    tts2 = getTravelTimes(dist, 0)
    tts = [tt for tt in tts if tt['phase_name'] == phase]
    tts2 = [tt for tt in tts2 if tt['phase_name'] == phase]
    if len(tts) == 0 or len(tts2) == 0:
        raise Exception('Taup does not return phase %s at event distance %s' %
                        (phase, dist))
    onset = event.origins[0].time + tts[0]['time']
    inc = tts2[0]['take-off angle']  # approximation
    return AttribDict({'dist':dist, 'back_azimuth':baz, 'inclination': inc,
                       'onset':onset})

def rf(method='dmt', *args, **kwargs):
    if method == 'dmt':
        rf_dmt(*args, **kwargs)
    elif method == 'client':
        rf_client(*args, **kwargs)

def rf_dmt(events='events_rf.xml', rftype='Ps', dist=(30, 90),
           **rf_kwargs):
    events = io.read_rfevents(events)
    print events
    for event in events:
        event_id = event.resource_id.getQuakeMLURI().split('/')[-1]
        output_file = os.path.join(conf.output_path, event_id, conf.rf_data)
        input_files = glob.glob(os.path.join(conf.dmt_path, event_id,
                                             conf.dmt_data))
        while len(input_files) > 0:
            files_tmp = input_files[0][:-1] + '?'
            for f in glob.glob(files_tmp):
                input_files.remove(f)
            st = read(files_tmp, headonly=True)
            io.read_sac_header(st)
            stats = event2stats(st[0].stats.latitude, st[0].stats.longitude,
                                event, phase=rftype[0], dist_range=dist)
            if not stats:
                continue
            st = read(files_tmp)
            st.merge()
            if len(st) != 3:
                import warnings
                warnings.warn('Need 3 component seismograms. More or less '
                              'than three components for files %s' % files_tmp)
                continue
            io.read_sac_header(st)
            rf_stream(st, stats, **rf_kwargs)
            io.write_sac_header(st, event)
            for tr in st:
                io.create_dir(output_file)
                tr.write(output_file.format(stats=tr.stats), 'SAC')

def rf_client(getwaveform, stations, events='events_rf.xml',
               request_window=(-50, 150), rftype='Ps', dist=(30, 90),
               **rf_kwargs):
    events = io.read_rfevents(events)
    print events
    if isinstance(stations, basestring):
        stations = io.read_stations(stations)
    for event in events:
        event_id = event.resource_id.getQuakeMLURI().split('/')[-1]
        output_file = os.path.join(conf.output_path, event_id, conf.rf_data)
        for station in stations:
            lat = stations[station].latitude
            lon = stations[station].longitude
            el = stations[station].elevation
            print lat, lon, event
            stats = event2stats(lat, lon, event, phase=rftype[0],
                                dist_range=dist)
            print stats
            if not stats:
                continue
            st = getwaveform(station, stats.onset + request_window[0],
                             stats.onset + request_window[1])
            print st
            st.merge()
            if len(st) != 3:
                import warnings
                warnings.warn('Need 3 component seismograms. More or less '
                              'than three components for event %s, station %s.'
                              % (event_id, station))
                continue
            stats.latitude = lat
            stats.longitude = lon
            stats.elevation = stations[station].elevation
            rf_stream(st, stats, **rf_kwargs)
            io.write_sac_header(st, event)
            io.create_dir(output_file)
            for tr in st:
                tr.write(output_file.format(stats=tr.stats), 'SAC')

def rf_stream(stream, stats, method='P', window=(-20, 100),
              downsample=None, filter=None,  #@ReservedAssignment
              rotate='ZNE->LQT', angles='stats', signoise=None,
              deconvolve='time', **deconvolve_kwargs):
    st = stream
    assert len(st) == 3
    if method != 'P':
        raise NotImplementedError
    for tr in st:
        tr.stats.update(stats)
    if filter:
        st.filter(filter)
    st.trim(stats.onset + window[0], stats.onset + window[1])
    if downsample and downsample <= st[0].stats.sampling_rate:
        st.decimate(int(st[0].stats.sampling_rate) // downsample)
    if angles == 'polar':
        raise NotImplementedError
    st.rotate(rotate)
    src_comp = rotate.split('->')[-1][0]
    if signoise:
        raise NotImplementedError
    if deconvolve:
        deconv(st, src_comp, method=deconvolve, **deconvolve_kwargs)
    # multiply -1 on Q and T component
    for tr in st:
        if tr.stats.channel[-1] != src_comp:
            tr.data = -tr.data
