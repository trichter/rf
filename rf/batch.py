"""
Functions for massive receiver function calculation. Untestet, in development.
"""
import argparse
import glob
import os
from obspy import read, readEvents
from rf.rfstream import rfstats, RFStream

try:
    import conf
    from rf import io
except ImportError:
    pass
#import warnings
#warnings.warn("Didn't find file conf.py. Using test configuration.")
#import rf.conf_test as conf






def rf_batch(method='dmt', *args, **kwargs):
    """
    TODO: doc rf_batch
    """
    if ' dist' not in kwargs:
        kwargs['dist'] = ((30, 90) if kwargs.get('method', 'P') == 'P' else
                          (60, 85))
    if method == 'dmt':
        rf_dmt(*args, **kwargs)
    elif method == 'client':
        rf_client(*args, **kwargs)


def rf_dmt(events=None, method='P', dist=None,
           **rf_kwargs):
    """
    TODO: doc rf_dmt
    """
    events = events or readEvents(conf.events)
    print events
    for event in events:
        event_id = event.resource_id.getQuakeMLURI().split('/')[-1]
        inputs = conf.data.format(eventid=event_id)
        inputs = glob.glob(os.path.join(conf.data_path, conf.data))
        while len(inputs) > 0:
            files_tmp = inputs[0][:-1] + '?'
            for f in glob.glob(files_tmp):
                inputs.remove(f)
            st = RFStream(read(files_tmp, headonly=True))
            st.read_sac_header()
            stats = rfstats(stats=st[0].stats, event=event, phase=method,
                            dist_range=dist)
            if not stats:
                continue
            st = RFStream(read(files_tmp))
            st.merge()
            if len(st) != 3:
                import warnings
                warnings.warn('Need 3 component seismograms. '
                              'Error for files %s' % files_tmp)
                continue
            for tr in st:
                tr.stats.update(stats)
            st.rf(**rf_kwargs)
            for tr in st:
                output = os.path.join(conf.output_path, conf.rf)
                output = output.format(eventid=event_id, stats=tr.stats)
                io._create_dir(output)
                tr.write(output, 'SAC')


def rf_client(getwaveform, stations=None, events=None,
              request_window=(-50, 150), method='P', dist=None, **rf_kwargs):
# S: -300 bis 300
    """
    TODO: doc rf_client
    """
    events = events or readEvents(conf.events)
    stations = stations or io.read_stations(conf.stations)
    for event in events:
        event_id = event.resource_id.getQuakeMLURI().split('/')[-1]
        for station in stations:
            stats = rfstats(station=stations[station], event=event,
                            phase=method, dist_range=dist)
            if not stats:
                continue
            st = getwaveform(station, stats.onset + request_window[0],
                             stats.onset + request_window[1])
            st = RFStream(stream=st)
            st.merge()
            if len(st) != 3:
                import warnings
                warnings.warn('Need 3 component seismograms. More or less '
                              'than three components for event %s, station %s.'
                              % (event_id, station))
                continue
            for tr in st:
                tr.stats.update(stats)
            st.rf(**rf_kwargs)
            st.write_sac_header()
            for tr in st:
                output = os.path.join(conf.output_path, conf.rf)
                output = output.format(eventid=event_id, stats=tr.stats)
                io._create_dir(output)
                tr.write(output, 'SAC')

def main():
    parser = argparse.ArgumentParser()
    print('Batch utility is not yet implemented.')

if __name__ == '__main__':
    main()

