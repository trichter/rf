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
from rf.rfstream import read_rf
from rf.util import _get_stations, iter_event_data, iter_event_metadata
try:
    from tqdm import tqdm
except ImportError:
    def tqdm():
        return None

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

PROFILE_FNAMES = {
    'Q': join('{root}', 'profile.QHD'),
    'SAC': join('{root}', 'profile_{box_pos}.SAC'),
    'H5': '{root}.h5'}


class _DummyDateTime(object):

    def __format__(self, *args, **kwargs):
        return '*'


class _DummyUTC(object):

    """Dummy UTCDateTime class returning '*' when formating"""

    def __init__(self):
        self.datetime = _DummyDateTime()


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


def _write(stream, root, format, type=None):
    if len(stream) == 0:
        return
    fname_pattern = (STACK_FNAMES if type == 'stack' else
                     PROFILE_FNAMES if type == 'profile' else
                     FNAMES)[format]
    fname = fname_pattern.format(root=root, **stream[0].stats)
    _create_dir(fname)
    if format == 'H5':
        stream.write(fname, format, mode='a')
    elif format == 'Q':
        stream.write(fname, format)
    elif format == 'SAC':
        for tr in stream:
            fname = fname_pattern.format(root=root, **tr.stats)
            tr.write(fname, format)


def _iter_event_processed_data(events, inventory, pin, format,
                               yield_traces=False, pbar=None):
    for meta in iter_event_metadata(events, inventory, pbar=pbar):
        meta['channel'] = '???'
        if 'event_time' not in meta and format != 'H5':
            meta['event_time'] = _DummyUTC()
        fname = FNAMES[format].format(root=pin, **meta)
        kwargs = {}
        if format == 'H5':
            meta.pop('channel')
            kwargs['readonly'] = meta
        try:
            stream = read_rf(fname, format, **kwargs)
        except:
            pass
        else:
            if yield_traces:
                for tr in stream:
                    yield tr
            else:
                yield stream


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
    try:
        client_module = import_module('obspy.clients.%s' % data)
    except ImportError:
        client_module = None
    if client_module:
        Client = getattr(client_module, 'Client')
        client = Client(**client_options)

        def get_waveforms(event=None, **args):
            return client.get_waveforms(**args)
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
            msg = 'channel %s: error while retrieving data: %s'
            log.debug(msg, seedid, ex)

    use_cache = client_module is not None or data == 'plugin'
    use_cache = use_cache and cache_waveforms
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


def run(args, conf=None, tutorial=False, get_waveforms=None, format='Q',
        **kwargs):
    command = args[0]
    # Create example configuration file and tutorial
    if command == 'create':
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
    # Load configuration
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
        # Populate kwargs with conf, but prefer kwargs
        conf.update(kwargs)
        kwargs = conf
        if 'moveout_phase' in kwargs:
            kwargs.setdefault('moveout', {})
            kwargs['moveout']['phase'] = kwargs.pop('moveout_phase')
    # Configure logging
    kw = {'loggingc': kwargs.pop('logging', None),
          'verbose': kwargs.pop('verbose', 0),
          'loglevel': kwargs.pop('loglevel', 3),
          'logfile': kwargs.pop('logfile', None)}
    configure_logging(**kw)
    # Read events and inventory
    try:
        if command in ('stack', 'plot'):
            events = None
        elif command != 'print' or args[1] == 'events':
            events = kwargs.pop('events')
            if (not isinstance(events, obspy.Catalog) or
                    not isinstance(events, list) or
                    (len(events) == 2 and isinstance(events[0], str))):
                if isinstance(events, (basestring, str)):
                    format_ = None
                else:
                    events, format_ = events
                events = obspy.read_events(events, format_)
                log.info('read %d events', len(events))
        if command != 'print' or args[1] == 'stations':
            inventory = kwargs.pop('inventory')
            if not isinstance(inventory, obspy.Inventory):
                if isinstance(inventory, (basestring, str)):
                    format_ = None
                else:
                    inventory, format_ = inventory
                inventory = obspy.read_inventory(inventory, format_)
                log.info('read inventory with %d stations',
                         len(_get_stations(inventory)))
    except (KeyboardInterrupt, SystemExit):
        raise
    except:
        log.exception('cannot read events or stations')
        return
    # Initialize get_waveforms
    if command == 'data':
        try:
            keys = ['client_options', 'plugin', 'cache_waveforms']
            tkwargs = {k: kwargs.pop(k, None) for k in keys}
            # Initialize get_waveforms
            if get_waveforms is None:
                data = kwargs.pop('data')
                get_waveforms = init_data(data, **tkwargs)
                log.info('init data from %s', data)
        except (KeyboardInterrupt, SystemExit):
            raise
        except:
            log.exception('cannot initalize data')
            return
    # Print command
    if command == 'print':
        if args[1] == 'events':
            print(events.__str__(True))
        elif args[1] == 'stations':
            print(inventory)
        else:
            raise NotImplementedError
        return
    # Dispatch positional arguments
    commands = []
    try:
        if command == 'data':
            pin = get_waveforms
            pout = args[-1]
            commands = args[:-1]
        elif command in ('calc', 'moveout'):
            pin, pout = args[-2:]
            commands = args[:-2]
        elif command in ('stack', 'plot', 'profile'):
            pin, pout = args[1:]
        elif command == 'convert':
            pin, pout, newformat = args[1:]
        else:
            raise
        for c in commands:
            assert c in ('data', 'calc', 'moveout')
    except:
        raise ParseError('Positional arguments not matching. Consult rf -h.')
    # Select appropriate iterator
    if command == 'data':
        options = kwargs.get('options', {})
        options.setdefault('phase', 'P')
        iter_ = iter_event_data(events, inventory, pin, pbar=tqdm(), **options)
    else:
        yt = command == 'profile'
        iter_ = _iter_event_processed_data(
            events, inventory, pin, format, pbar=tqdm(), yield_traces=yt)
    # Run all subcommands
    if command == 'convert':
        for stream in iter_:
            _write(stream, pout, newformat)
    elif command == 'plot':
        for stream in iter_:
            channels = set(tr.stats.channel for tr in stream)
            for ch in channels:
                st2 = stream.select(channel=ch)
                fname = PLOT_FNAMES.format(root=pout, **st2[0].stats)
                _create_dir(fname)
                st2.sort(['back_azimuth'])
                st2.plot_rf(fname, **kwargs.get('plot', {}))
    elif command == 'stack':
        for stream in iter_:
            stack = stream.stack()
            _write(stack, pout, format, type='stack')
    elif command == 'profile':
        from rf.profile import get_profile_boxes, get_profile
        boxes = get_profile_boxes(**kwargs.get('boxes', {}))
        profile = get_profile(iter_, boxes, **kwargs.get('profile', {}))
        _write(profile, pout, format, type='profile')
    else:
        for stream in iter_:
            for command in commands:
                if command == 'data':
                    pass
                elif command == 'calc':
                    stream.rf(**kwargs.get('rf', {}))
                elif command == 'moveout':
                    stream.moveout(**kwargs.get('moveout', {}))
                else:
                    raise NotImplementedError
            _write(stream, pout, format)


def run_cli(args=None):
    from rf import __version__
    p = argparse.ArgumentParser(description=__doc__)
    version = '%(prog)s ' + __version__
    p.add_argument('--version', action='version', version=version)
    msg = 'Configuration file to load (default: conf.json)'
    p.add_argument('-c', '--conf', default='conf.json', help=msg)
    msg = 'Set chattiness on command line. Up to 3 -v flags are possible'
    p.add_argument('-v', '--verbose', help=msg, action='count',
                   default=SUPPRESS)
    p.add_argument('args', nargs='+')

#
#
#    sub = p.add_subparsers(title='subcommands', dest='subcommand')
#    msg = 'create config file in current directory'
#    p_create = sub.add_parser('create', help=msg)
#    msg = 'create some example files, too'
    p.add_argument('-t', '--tutorial', help=msg, action='store_true')
#
#    msg = 'calculate receiver functions'
#    p_calc = sub.add_parser('calc', help=msg)
#    msg = 'print information about events or stations'
#    p_print = sub.add_parser('print', help=msg)
#    msg = 'print information about subject'
#    p_print.add_argument('subject', help=msg, choices=('stations', 'events'))
#    msg = 'perform move out correction'
#    p_mout = sub.add_parser('moveout', help=msg)
#    msg = 'convert files to different format'
#    p_conv = sub.add_parser('convert', help=msg)
#    msg = 'plot receiver functions'
#    p_plot = sub.add_parser('plot', help=msg)
#    msg = 'stack receiver functions'
#    p_stack = sub.add_parser('stack', help=msg)
#
#    for pp in [p_conv, p_plot, p_stack]:
#        pp.add_argument('root', help='directory or basename of file')
#    msg = 'convert to other format (supported: Q, SAC or H5)'
#    p_conv.add_argument('newformat', help=msg)

    msg = ('Use these flags to overwrite values in the config file. '
           'See the example configuration file for a description of '
           'these options')

    g2 = p.add_argument_group('optional config arguments', description=msg)
    features_str = ('events', 'inventory', 'data', 'phase',
                    'moveout_phase')
    for f in features_str:
        g2.add_argument('--' + f, default=SUPPRESS)

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
