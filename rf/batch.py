"""
Functions for massive receiver function calculation. Untestet, in development.
"""
import argparse
import glob
import os
import shutil
from obspy import read, readEvents
from rf.rfstream import rfstats, RFStream

##### start ex rf.io #####
import pickle
from obspy.core.event import (Catalog, Event, CreationInfo, EventDescription,
                              Origin, Magnitude)
from obspy.core.util import AttribDict

#def _manipulate_conf(conf):
#    conf.rf_in = conf.rf
#    conf.mout_in = conf.mout
#    for key, val in CONF_ABBREVS.items():
#        conf.data = conf.data.replace(key, '*')
#        conf.rf_in = conf.rf_in.replace(key, '*')  # @UndefinedVariable
#        conf.mout_in = conf.mout_in.replace(key, '*')  # @UndefinedVariable
#        conf.rf = conf.rf.replace(key, val)
#        conf.mout = conf.mout.replace(key, val)
#        conf.mean = conf.mean.replace(key, val)


def _create_dir(filename):
    """
    Create directory of filname if it does not exist.
    """
    head = os.path.dirname(filename)
    if head != '' and not os.path.isdir(head):
        os.makedirs(head)


def _read_stations(fname):
    """
    Read station positions from whitespace delimited file

    Example file:
    # station  lat  lon  elev
    STN  10.0  -50.0  160
    """
    ret = AttribDict()
    with open(fname) as f:
        for line in f.readlines():
            if not line[0].startswith('#'):
                vals = line.split()
                ret[vals[0]] = AttribDict()
                ret[vals[0]].latitude = float(vals[1])
                ret[vals[0]].longitude = float(vals[2])
                ret[vals[0]].elevation = float(vals[3])
    return ret


def _write_stations(fname, dic, comment='# station  lat  lon  elev\n'):
    """
    Write dictionary of station positions to file
    """
    with open(fname, "w") as f:
        f.write(comment)
        for key, val in sorted(dic.items()):
            f.write('%s  %s  %s  %s\n' % (key, val['latitude'],
                                          val['longitude'], val['elevation']))


#def _convert_dmteventfile():
#    eventsfile1 = os.path.join(conf.data_path, 'EVENT', 'event_list')
#    eventsfile2 = os.path.join(conf.data_path, 'EVENT', 'events.xml')
#    with open(eventsfile1) as f:
#        events1 = pickle.load(f)
#    events2 = Catalog()
#    for ev in events1:
#        orkw = {'time': ev['datetime'],
#                'latitude': ev['latitude'],
#                'longitude': ev['longitude'],
#                'depth': ev['depth']}
#        magkw = {'mag': ev['magnitude'],
#                 'magnitude_type': ev['magnitude_type']}
#        evdesargs = (ev['flynn_region'], 'Flinn-Engdahl region')
#        evkw = {'resource_id': ev['event_id'],
#                'event_type': 'earthquake',
#                'creation_info': CreationInfo(author=ev['author']),
#                'event_descriptions': [EventDescription(*evdesargs)],
#                'origins': [Origin(**orkw)],
#                'magnitudes': [Magnitude(**magkw)]}
#        events2.append(Event(**evkw))
#    events2.write(eventsfile2, 'QUAKEML')
#
#def _create_rfeventsfile(events='events.xml',
#                        eventsfile='events_rf.xml', filters=None):
#    if isinstance(events, basestring):
#        if not os.path.exists(events):
#            events = os.path.join(conf.data_path, 'EVENT', 'events.xml')
#        events = readEvents(events)
#    if filters:
#        events.filter(*filters)
#    eventsfile = os.path.join(conf.output_path, 'EVENT', eventsfile)
#    _create_dir(eventsfile)
#    events.write(eventsfile, 'QUAKEML')

##### end ex rf.io #####


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


def rf_dmt(data_path, rf, events=None, phase='P', dist=None,
           **rf_kwargs):
    """
    TODO: doc rf_dmt
    """
    events = readEvents(events)
    print events
    for event in events:
        event_id = event.resource_id.getQuakeMLURI().split('/')[-1]
        inputs = data_path.format(eventid=event_id)
        inputs = glob.glob(data_path)
        while len(inputs) > 0:
            files_tmp = inputs[0][:-1] + '?'
            for f in glob.glob(files_tmp):
                inputs.remove(f)
            st = RFStream(read(files_tmp, headonly=True))
            st.read_sac_header()
            stats = rfstats(stats=st[0].stats, event=event, phase=phase,
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
            st.rf(method=phase[0], **rf_kwargs)
            for tr in st:
                output = rf.format(eventid=event_id, stats=tr.stats)
                _create_dir(output)
                tr.write(output, 'SAC')


def rf_client(get_waveform, rf, stations=None, events=None,
              request_window=(-50, 150), phase='P', dist=None,
              **rf_kwargs):
# S: -300 bis 300
    """
    TODO: doc rf_client
    """
    events = readEvents(events)
    stations = _read_stations(stations)
    for event in events:
        event_id = event.resource_id.getQuakeMLURI().split('/')[-1]
        for station in stations:
            stats = rfstats(station=stations[station], event=event,
                            phase=phase, dist_range=dist)
            if not stats:
                continue
            st = get_waveform(station, stats.onset + request_window[0],
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
            st.rf(method=phase[0], **rf_kwargs)
            st.write_sac_header()
            for tr in st:
                output = rf.format(eventid=event_id, stats=tr.stats)
                _create_dir(output)
                tr.write(output, 'SAC')


def _filter_config(config):
    method = config['method']
    deconvolve = config['deconvolve']
    settings = ['method', 'phase', 'events', 'stations', 'rf', 'mout', 'mean',
               'format', 'filter', 'window', 'downsample', 'rotate',
               'deconvolve']
    options = ['winsrc']
    if method == 'dmt':
        settings.append('data_path')
    else:
        settings.extend(['get_waveform', 'request_window'])
    if deconvolve == 'freq':
        settings.extend(['water', 'gauss'])
        options.append('tshift')
    else:
        settings.append('spiking')
        options.extend(['winrsp', 'winrf'])
    for key in config:
        if key not in settings and key not in options:
            config.pop(key)

def main():
    parser = argparse.ArgumentParser(description=
                                     'rf batch command line utility')
    subparsers = parser.add_subparsers(title='subcommands', dest='subcommand')
    parser_init = subparsers.add_parser('init', help='initialize new project')
    parser_init.add_argument('path', help='path for new rf project')
    parser_calc = subparsers.add_parser('calc', help='calculate receiver '
                                                     'functions')
    parser_calc.add_argument('-c', '--conf', default='conf.py',
                             help='config file name')
    args = parser.parse_args()
    print('Batch utility is not yet implemented.')
    return
    if args.subcommand == 'init':
        path = args.path
        if os.path.exists(path):
            print(('Directory %s exists already. To create a new rf project '
                   'please delete the directory beforehand.') % path)
            return
        os.mkdir(path)
        src = os.path.join(os.path.dirname(__file__), 'conf_test.py')
        dst = os.path.join(path, 'conf.py')
        shutil.copyfile(src, dst)
        print('New rf project initialized in directory %s' % path)
        return
    conf = {}
    execfile(args.conf, conf)
    try:
        conf['window'] = conf['window' + conf['phase'][0].upper()]
    except KeyError:
        pass
    _filter_config(conf)
    if conf.pop('method') == 'dmt':
        rf_dmt(**conf)
    else:
        rf_client(**conf)


if __name__ == '__main__':
    main()

