"""
rf: receiver function calculation - batch command line utility
"""

import argparse
from argparse import SUPPRESS
from copy import deepcopy
from importlib import import_module
import json
import logging
import logging.config
import os
from os.path import join
from pkg_resources import resource_filename
import shutil
import sys

import obspy
from rf.rfstream import read_rf, rfstats, RFStream, set_index
try:
    from progressbar import ProgressBar
except ImportError:
    ProgressBar = None

try:
    import joblib
except ImportError:
    joblib = None

IS_PY3 = sys.version_info.major == 3

log = logging.getLogger('rf.batch')
log.addHandler(logging.NullHandler())

LOGLEVELS = {0: 'CRITICAL', 1: 'WARNING', 2: 'INFO', 3: 'DEBUG'}

LOGGING_DEFAULT_CONFIG = {
    'version': 1,
    'disable_existing_loggers': False,
    'capture_warnings': True,
    'formatters': {
        'file': {
            'format': ('%(asctime)s %(module)-6s%(process)-6d%(levelname)-8s'
                       '%(message)s')
        },
        'console': {
            'format': '%(levelname)-8s%(message)s'
        },
    },
    'handlers': {
        'console': {
            'class': 'logging.StreamHandler',
            'formatter': 'console',
            'level': None,
        },
        'file': {
            'class': 'logging.FileHandler',
            'formatter': 'file',
            'level': None,
            'filename': None,
        },
    },
    'loggers': {
        'rf': {
            'handlers': ['console', 'file'],
            'level': 'DEBUG',
            'propagate': False,
        },
        'py.warnings': {
            'handlers': ['console', 'file'],
            'level': 'DEBUG',
            'propagate': False,
        }

    }
}

_TF = '.datetime:%Y-%m-%dT%H:%M:%S'

FNAMES = {
    'Q': join('{root}', '{network}.{station}.{location}',
              '{network}.{station}.{location}_{event_time%s}.QHD' % _TF),
    'SAC': join('{root}', '{network}.{station}.{location}',
                '{network}.{station}.{location}.{channel}_'
                '{event_time%s}.SAC' % _TF),
    'H5': '{root}.h5'}
PLOT_FNAMES = join('{root}', '{network}.{station}.{location}.{channel}.pdf')
STACK_FNAMES = {
    'Q': join('{root}', '{network}.{station}.{location}.QHD'),
    'SAC': join('{root}', '{network}.{station}.{location}.{channel}.SAC'),
    'H5': '{root}.h5'}


def _no_pbar():
    """Turn off progressbar"""
    global ProgressBar
    ProgressBar = None


class _DummyDateTime(object):

    def __format__(self, *args, **kwargs):
        return '*'


class _DummyUTC(object):

    """Dummy UTCDateTime class returning '*' when formating"""

    def __init__(self):
        self.datetime = _DummyDateTime()


def _check_path(path):
    if os.path.exists(path):
        print(('Directory %s exists already. You have to delete it '
               'beforehand to perform this action.') % path)
        sys.exit()


def _create_dir(filename):
    """
    Create directory of filname if it does not exist.
    """
    head = os.path.dirname(filename)
    if head != '' and not os.path.isdir(head):
        os.makedirs(head)


def _fname(fname, **kwargs):
    fname = fname.format(**kwargs)
    _create_dir(fname)
    return fname


def _write(stream, path, root, format, stack=False):
    fnames = STACK_FNAMES if stack else FNAMES
    fname = join(path, fnames[format])
    if format == 'H5':
        set_index('rf_stack' if stack else 'rf')
        stream.write(_fname(fname, root=root), format, mode='a')
    elif format == 'Q':
        stream.write(_fname(fname, root=root, **stream[0].stats), format)
    elif format == 'SAC':
        for tr in stream:
            tr.write(_fname(fname, root=root, **tr.stats), format)


def _read(stats, path, root, format):
    fname = os.path.join(path, FNAMES[format])
    if format == 'H5':
        stats.pop('channel')
        if isinstance(stats['event_time'], _DummyUTC):
            stats.pop('event_time')
    fname = fname.format(root=root, **stats)
    kwargs = {}
    if format == 'H5':
        set_index()
        kwargs['readonly'] = stats
    try:
        return read_rf(fname, format, **kwargs)
    except:
        pass


def _iter(*args, **kwargs):
    """Return _Iter iterable with optional progressbar support"""
    if ProgressBar:
        return ProgressBar()(_Iter(*args, **kwargs))
    else:
        return _Iter(*args, **kwargs)


class _Iter(object):

    """Iterable wrapping _generator with length"""

    def __init__(self, *args, **kwargs):
        self.gen = _generator(*args, **kwargs)
        self.length = next(self.gen)

    def __len__(self):
        return self.length

    def __iter__(self):
        return self.gen


def _generator(events, inventory, rf=False):
    """Generator yielding length at first and then station/event information"""
    channels = inventory.get_contents()['channels']
    stations = list(set(ch.rsplit('.', 1)[0] for ch in channels))
    one_channel = {ch.rsplit('.', 1)[0]: ch for ch in channels}
    if events is not None:
        yield len(stations) * len(events)
        for event in events:
            for station in stations:
                seed_id = one_channel[station][:-1] + '?'
                net, sta, loc, cha = seed_id.split('.')
                stats = {'network': net, 'station': sta, 'location': loc,
                         'channel': cha}
                if rf:
                    stats['event'] = event
                    #stats['seed_id'] = seed_id
                    coords = inventory.get_coordinates(one_channel[station])
                    yield stats, event, coords
                else:
                    stats['event_time'] = event.preferred_origin()['time']
                    yield stats
    else:
        yield len(stations)
        for station in stations:
            net, sta, loc, cha = one_channel[station].split('.')
            stats = {'network': net, 'station': sta, 'location': loc,
                     'channel': cha[:-1] + '?',
                     'event_time': _DummyUTC()}
            yield stats


def run_rf(events, inventory, get_waveforms, path, format='H5',
           request_window=None, phase='P', dist_range=None,
           **rf_kwargs):
    root = phase + 'rf'
    _check_path(join(path, root))
    method = phase[-1].upper()
    if dist_range is None:
        dist_range = (30, 90) if method == 'P' else (60, 85)
    if request_window is None:
        request_window = (-50, 150) if method == 'P' else (-100, 50)
    for kwargs, event, coords in _iter(events, inventory, rf=True):
        stats = rfstats(station=coords, event=event,
                        phase=phase, dist_range=dist_range)
        if not stats:
            continue
        kwargs.update({'starttime': stats.onset + request_window[0],
                       'endtime': stats.onset + request_window[1]})
        stream = get_waveforms(**kwargs)
        if stream is None:
            continue
        stream = RFStream(stream, warn=False)
        stream.merge()
        if len(stream) != 3:
            import warnings
            warnings.warn('Need 3 component seismograms. More or less '
                          'than three components for event %s, station %s.'
                          % (stats.event_id, kwargs['seed_id']))
            continue
        for tr in stream:
            tr.stats.update(stats)
        stream.rf(method=method, **rf_kwargs)
        if len(stream) != 3:
            continue
        _write(stream, path, root, format)


def run_moveout(moveout, events, inventory, format,
                kwargs_moveout, path, root=None, phase='P', **kwargs):
    if root is None:
        root = phase[-1].upper() + 'rf'
    for stats in _iter(events, inventory):
        stream = _read(stats, path, root, format)
        if stream is None:
            continue
        stream.moveout(moveout, **kwargs_moveout)
        if len(stream) > 0:
            _write(stream, path, '%s_%s' % (root, moveout), format)


def run_convert(events, inventory, root, newformat, format, path='.', **kwargs):
    for stats in _iter(events, inventory):
        stream = _read(stats, path, root, format)
        if stream is None:
            continue
        _write(stream, path, join(newformat, root), newformat)


def run_plot(events, inventory, root, format, kwargs_plot, path='.',
             **kwargs):
    for stats in _iter(None, inventory):
        stream = _read(stats, path, root, format)
        if stream is None:
            continue
        channels = set(tr.stats.channel for tr in stream)
        for ch in channels:
            st2 = stream.select(channel=ch)
            stats['channel'] = ch
            fname = PLOT_FNAMES.format(root='plot_%s' % root, **stats)
            fname = join(path, fname)
            _create_dir(fname)
            st2.sort(['back_azimuth'])
            st2.plot_rf(fname, **kwargs_plot)


def run_stack(inventory, path, root, format, **kwargs):
    for stats in _iter(None, inventory):
        stream = _read(stats, path, root, format)
        stream.stack()
        _write(stream, path, 'stack_%s' % root, format, stack=True)

CONFIG = ['events', 'inventory', 'get_waveforms', 'request_window',
          'dist_range', 'path', 'format',
          'phase', 'filter', 'window', 'downsample', 'rotate',
          'deconvolve', 'source_component', 'winsrc']
CONFIG_FREQ = ['water', 'gauss', 'tshift']
CONFIG_TIME = ['winrsp', 'winrf', 'spiking']

def _slice_config(conf):
    if conf['deconvolve'] == 'freq':
        CONFIG.extend(CONFIG_FREQ)
    else:
        CONFIG.extend(CONFIG_TIME)
    return {k: conf[k] for k in CONFIG if k in conf}


def load_func(modulename, funcname):
    """Load and return function from Python module"""
    sys.path.append(os.path.curdir)
    module = import_module(modulename)
    sys.path.pop(-1)
    func = getattr(module, funcname)
    return func


def init_data(data, client_options=None, plugin=None, cache_waveforms=False):
    """Return appropriate get_waveforms function

    See example configuration file for a description of the options"""
    if client_options is None:
        client_options = {}
    is_webservice = data in ('arclink', 'fdsn', 'seishub')
    if is_webservice:
        webservice_module = import_module('obspy.%s' % data)
        Client = getattr(webservice_module, 'Client')
        client = Client(**client_options)
        if data == 'fdsn':
            get_waveforms_orig = client.get_waveforms
        else:
            get_waveforms_orig = client.getWaveform

        def get_waveforms(event=None, **args):
            return get_waveforms_orig(**args)
    elif data == 'plugin':
        modulename, funcname = plugin.split(':')
        get_waveforms = load_func(modulename.strip(), funcname.strip())
    else:
        from obspy import read
        stream = read(data)

        def get_waveforms(network, station, location, channel,
                          starttime, endtime, event=None):
            st = stream.select(network=network, station=station,
                               location=location, channel=channel)
            st = st.slice(starttime, endtime)
            return st

    def wrapper(**kwargs):
        try:
            return get_waveforms(**kwargs)
        except (KeyboardInterrupt, SystemExit):
            raise
        except Exception as ex:
            seedid = '.'.join((kwargs['network'], kwargs['station'],
                               kwargs['location'], kwargs['channel']))
            msg = 'channel %s: error while retireving data: %s'
            log.debug(msg, seedid, ex)

    use_cache = cache_waveforms and (is_webservice or data == 'plugin')
    if use_cache and joblib:
        log.info('use waveform cache in %s', cache_waveforms)
        memory = joblib.Memory(cachedir=cache_waveforms, verbose=0)
        return memory.cache(wrapper)
    elif use_cache:
        log.warning('install joblib to use cache_waveforms option')
    return wrapper


class ConfigJSONDecoder(json.JSONDecoder):

    """Strip lines from comments"""

    def decode(self, s):
        s = '\n'.join(l.split('#', 1)[0] for l in s.split('\n'))
        return super(ConfigJSONDecoder, self).decode(s)

def configure_logging(loggingc, verbose=0, loglevel=3, logfile=None):
    if loggingc is None:
        loggingc = deepcopy(LOGGING_DEFAULT_CONFIG)
        if verbose > 3:
            verbose = 3
        loggingc['handlers']['console']['level'] = LOGLEVELS[verbose]
        if logfile is None or loglevel == 0:
            del loggingc['handlers']['file']
            loggingc['loggers']['qopen']['handlers'] = ['console']
            loggingc['loggers']['py.warnings']['handlers'] = ['console']
        else:
            loggingc['handlers']['file']['level'] = LOGLEVELS[loglevel]
            loggingc['handlers']['file']['filename'] = logfile
    logging.config.dictConfig(loggingc)
    logging.captureWarnings(loggingc.get('capture_warnings', False))

class ParseError(Exception):
    pass

def get_station(seedid):
    """Station name from seed id"""
    st = seedid.rsplit('.', 2)[0]
    if st.startswith('.'):
        st = st[1:]
    return st

#def get_eventid(event):
#    """Event id from event"""
#    return str(event.resource_id).split('/')[-1]
#
#
#def get_pair(tr):
#    """Station and event id from trace"""
#    return (tr.stats.eventid, get_station(tr.id))


def run(subcommand, conf=None, tutorial=False, get_waveforms=None,
        subject='events', **args):
    if subcommand == 'create':
        if conf is None:
            conf = 'conf.json'
        srcs = ['conf.json']
        dest_dir = os.path.dirname(conf)
        dests = [conf]
        if tutorial:
            example_files = ['example_events.xml', 'example_inventory.xml',
                             'example_data.mseed']
            srcs.extend(example_files)
            for src in example_files:
                dests.append(os.path.join(dest_dir, src))
        for src, dest in zip(srcs, dests):
            src = resource_filename('rf', 'example/%s' % src)
            shutil.copyfile(src, dest)
        return
    if conf in ('None', 'none', 'null', ''):
        conf = None
    if conf:
        try:
            with open(conf) as f:
                conf = json.load(f, cls=ConfigJSONDecoder)
        except ValueError as ex:
            print('Error while parsing the configuration: %s' % ex)
            return
        except IOError as ex:
            print(ex)
            return
        # Populate args with conf, but prefer args
        conf.update(args)
        args = conf
    # Configure logging
    kw = {'loggingc': args.pop('logging', None),
          'verbose': args.pop('verbose', 0),
          'loglevel': args.pop('loglevel', 3),
          'logfile': args.pop('logfile', None)}
    configure_logging(**kw)
    # Read events and inventory
    try:
        if subcommand != 'print' or subject == 'events':
            events = args.pop('events')
            if not isinstance(events, (list, obspy.core.event.Catalog)):
                events = obspy.readEvents(events)
                log.info('read %d events', len(events))
            args['events'] = events
        if subcommand != 'print' or subject == 'stations':
            inventory = args.pop('inventory')
            if not isinstance(inventory, obspy.station.Inventory):
                inventory = obspy.read_inventory(inventory)
                channels = inventory.get_contents()['channels']
                stations = list(set(get_station(ch) for ch in channels))
                log.info('read inventory with %d stations', len(stations))
            args['inventory'] = inventory
    except (KeyboardInterrupt, SystemExit):
        raise
    except:
        log.exception('cannot read events or stations')
        return
    # Run subcommands
    if subcommand == 'print':
        if subject == 'events':
            print(events.__str__(True))
        else:
            print(inventory)
    elif subcommand == 'calc':
        try:
            # Initialize get_waveforms
            keys = ['client_options', 'plugin', 'cache_waveforms']
            tkwargs = {k: args.pop(k, None) for k in keys}
            if get_waveforms is None:
                data = args.pop('data')
                get_waveforms = init_data(data, **tkwargs)
                log.info('init data from %s', data)
        except (KeyboardInterrupt, SystemExit):
            raise
        except:
            log.exception('cannot initalize data')
            return
        args['get_waveforms'] = get_waveforms
        phase = args.get('phase', 'P')
        method = phase[-1].upper()
        options_by_phase=args.pop('options_by_phase')
        if options_by_phase:
            args.update(options_by_phase[method])
        args = _slice_config(args)
        run_rf(**args)
    elif subcommand == 'moveout':
        run_moveout(**args)
    elif subcommand == 'convert':
        run_convert(**args)
    elif subcommand == 'plot':
        run_plot(**args)
    elif subcommand == 'stack':
        run_stack(**args)
    else:
        print('Subcommand not availlable')

def run_cli(args=None):
    p = argparse.ArgumentParser(description=__doc__)
    msg = 'Configuration file to load (default: conf.json)'
    p.add_argument('-c', '--conf', default='conf.json', help=msg)
    msg = 'Set chattiness on command line. Up to 3 -v flags are possible'
    p.add_argument('-v', '--verbose', help=msg, action='count',
                   default=SUPPRESS)
    sub = p.add_subparsers(title='subcommands', dest='subcommand')
    msg = 'create config file in current directory'
    p_create = sub.add_parser('create', help=msg)
    msg = 'create some example files, too'
    p_create.add_argument('-t', '--tutorial', help=msg, action='store_true')

    msg = 'calculate receiver functions'
    p_calc = sub.add_parser('calc', help=msg)
    msg = 'print information about events or stations'
    p_print = sub.add_parser('print', help=msg)
    msg = 'print information about subject'
    p_print.add_argument('subject', help=msg, choices=('stations', 'events'))
    msg = 'perform move out correction'
    p_mout = sub.add_parser('moveout', help=msg)
    msg = 'convert files to different format'
    p_conv = sub.add_parser('convert', help=msg)
    msg = 'plot receiver functions'
    p_plot = sub.add_parser('plot', help=msg)
    msg = 'stack receiver functions'
    p_stack = sub.add_parser('stack', help=msg)

    for pp in [p_conv, p_plot, p_stack]:
        pp.add_argument('root', help='directory or basename of file')
    msg = 'convert to other format (supported: Q, SAC or H5)'
    p_conv.add_argument('newformat', help=msg)

    msg = ('Use these flags to overwrite values in the config file. '
           'See the example configuration file for a description of '
           'these options')

    g1 = p.add_argument_group('optional rf arguments', description=msg)
    features_str = ('events', 'inventory', 'data', 'method', 'phase',
                    'moveout', 'format', 'path')
    for f in features_str:
        g1.add_argument('--' + f, default=SUPPRESS)

#    features_bool = ()
#    for f in features_bool:
#        g1.add_argument('--' + f.replace('_', '-'), dest=f,
#                        action='store_true', default=SUPPRESS)
#        g1.add_argument('--no-' + f.replace('_', '-'), dest=f,
#                        action='store_false', default=SUPPRESS)

    # Get command line arguments and start run
    args = vars(p.parse_args(args))
    try:
        run(**args)
    except ParseError as ex:
        p.error(ex)
