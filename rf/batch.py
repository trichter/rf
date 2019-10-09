# Copyright 2013-2019 Tom Eulenfeld, MIT license
"""
rf: receiver function calculation - batch command line utility
"""

import argparse
from argparse import SUPPRESS
from importlib import import_module
import json
import os
from os.path import join
from pkg_resources import resource_filename
import shutil
import sys

import numpy as np
import obspy
from rf.rfstream import read_rf
from rf.util import iter_event_data, iter_event_metadata

try:
    from tqdm import tqdm
except ImportError:
    def tqdm():
        return None

try:
    basestring
except NameError:
    basestring = str
IS_PY3 = sys.version_info.major == 3

_TF = '.datetime:%Y-%m-%dT%H:%M:%S'

FNAMES = {
    'Q': join('{root}', '{network}.{station}.{location}',
              '{network}.{station}.{location}_{event_time%s}.QHD' % _TF),
    'SAC': join('{root}', '{network}.{station}.{location}',
                '{network}.{station}.{location}.{channel}_'
                '{event_time%s}.SAC' % _TF),
    'H5': '{root}.h5'}
STACK_FNAMES = {
    'Q': join('{root}', '{network}.{station}.{location}.QHD'),
    'SAC': join('{root}', '{network}.{station}.{location}.{channel}.SAC'),
    'H5': '{root}.h5'}
PROFILE_FNAMES = {
    'Q': '{root}.QHD',
    'SAC': join('{root}', 'profile_{channel[2]}_{box_pos}.SAC'),
    'H5': '{root}.h5'}
PLOT_FNAMES = join('{root}', '{network}.{station}.{location}.{channel}.pdf')
PLOT_PROFILE_FNAMES = join('{root}', 'profile_{channel[2]}.pdf')


class _DummyDateTime(object):

    def __format__(self, *args, **kwargs):
        return '*'


class _DummyUTC(object):

    """Dummy UTCDateTime class returning '*' when formatting."""

    def __init__(self):
        self.datetime = _DummyDateTime()


def _create_dir(filename):
    """
    Create directory of filname if it does not exist.
    """
    head = os.path.dirname(filename)
    if head != '' and not os.path.isdir(head):
        os.makedirs(head)


def write(stream, root, format, type=None):
    """Write stream to one or more files depending on format."""
    format = format.upper()
    if len(stream) == 0:
        return
    fname_pattern = (STACK_FNAMES if type == 'stack' else
                     PROFILE_FNAMES if type == 'profile' else
                     FNAMES)[format]
    fname = fname_pattern.format(root=root, **stream[0].stats)
    _create_dir(fname)
    if format == 'H5':
        stream.write(fname, format, mode='a', ignore=('mseed',))
    elif format == 'Q':
        stream.write(fname, format)
    elif format == 'SAC':
        for tr in stream:
            fname = fname_pattern.format(root=root, **tr.stats)
            tr.write(fname, format)


def iter_event_processed_data(events, inventory, pin, format,
                              yield_traces=False, pbar=None):
    """Iterator yielding streams or traces which are read from disc."""
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
        except Exception:
            pass
        else:
            if yield_traces:
                for tr in stream:
                    yield tr
            else:
                yield stream


def _iter_profile(pin, format):
    fname = PROFILE_FNAMES[format].format(root=pin, box_pos='*', channel='???')
    yield read_rf(fname)


def load_func(modulename, funcname):
    """Load and return function from Python module."""
    sys.path.append(os.path.curdir)
    module = import_module(modulename)
    sys.path.pop(-1)
    func = getattr(module, funcname)
    return func


def init_data(data, client_options=None, plugin=None):
    """Return appropriate get_waveforms function.

    See example configuration file for a description of the options."""
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
            print(msg % (seedid, ex))

    return wrapper


class ConfigJSONDecoder(json.JSONDecoder):
    """Strip lines from comments."""

    def decode(self, s):
        s = '\n'.join(l.split('#', 1)[0] for l in s.split('\n'))
        return super(ConfigJSONDecoder, self).decode(s)


class ParseError(Exception):
    pass


def run(command, conf=None, tutorial=False, **kw):
    """Create example configuration file and tutorial or load config.

    After that call `run_commands`.
    """
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
    if conf and (command != 'print' or
                 kw.get('objects', [''])[0] in ('stations', 'events')):
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
        conf.update(kw)
        kw = conf
    run_commands(command, **kw)


DICT_OPTIONS = ['client_options', 'options', 'rf', 'moveout',
                'boxbins', 'boxes', 'profile', 'plot', 'plot_profile']


def run_commands(command, commands=(), events=None, inventory=None,
                 objects=None, get_waveforms=None, data=None, plugin=None,
                 phase=None, moveout_phase=None,
                 path_in=None, path_out=None, format='Q',
                 newformat=None, **kw):
    """Load files, apply commands and write result files."""
    for opt in kw:
        if opt not in DICT_OPTIONS:
            raise ParseError('Unknown config option: %s' % opt)
    for opt in DICT_OPTIONS:
        default = None if opt == 'boxbins' else {}
        d = kw.setdefault(opt, default)
        if isinstance(d, basestring):
            kw[opt] = json.loads(d)
    if phase is not None:
        kw['options']['phase'] = phase
    if moveout_phase is not None:
        kw['moveout']['phase'] = moveout_phase
    if kw['boxbins'] is not None:
        kw['boxes']['bins'] = np.linspace(*kw['boxbins'])

    try:
        if command == 'calc':
            assert len(commands) < 3
            if len(commands) == 2:
                assert commands[0] != commands[1]
        elif command == 'calc':
            assert len(commands) < 2
    except Exception:
        raise ParseError('calc or moveout command given more than once')

    # Read events and inventory
    try:
        if command in ('stack', 'plot'):
            events = None
        elif command != 'print' or objects[0] == 'events':
            if (not isinstance(events, obspy.Catalog) or
                    not isinstance(events, list) or
                    (len(events) == 2 and isinstance(events[0], basestring))):
                if isinstance(events, basestring):
                    format_ = None
                else:
                    events, format_ = events
                events = obspy.read_events(events, format_)
        if command != 'print' or objects[0] == 'stations':
            if not isinstance(inventory, obspy.Inventory):
                if isinstance(inventory, basestring):
                    format_ = None
                else:
                    inventory, format_ = inventory
                inventory = obspy.read_inventory(inventory, format_)
    except Exception:
        print('cannot read events or stations')
        return
    # Initialize get_waveforms
    if command == 'data':
        try:
            # Initialize get_waveforms
            if get_waveforms is None:
                get_waveforms = init_data(
                    data, client_options=kw['client_options'], plugin=plugin)
        except Exception:
            print('cannot initalize data')
            return
    # Print command
    if command == 'print':
        if objects[0] == 'events':
            print(events.__str__(True))
        elif objects[0] == 'stations':
            print(inventory)
        else:
            from rf.rfstream import RFStream
            stream = sum((read_rf(fname) for fname in objects), RFStream())
            print(stream.__str__(True))
        return
    # Select appropriate iterator
    if command == 'data':
        iter_ = iter_event_data(events, inventory, get_waveforms, pbar=tqdm(),
                                **kw['options'])
    elif command == 'plot-profile':
        iter_ = _iter_profile(path_in, format)
    else:
        yt = command == 'profile'
        iter_ = iter_event_processed_data(
            events, inventory, path_in, format, pbar=tqdm(), yield_traces=yt)
    # Run all commands
    if command == 'convert':
        for stream in iter_:
            write(stream, path_out, newformat)
    elif command == 'plot':
        for stream in iter_:
            channels = set(tr.stats.channel for tr in stream)
            for ch in channels:
                st2 = stream.select(channel=ch)
                fname = PLOT_FNAMES.format(root=path_out, **st2[0].stats)
                _create_dir(fname)
                st2.sort(['back_azimuth'])
                st2.plot_rf(fname, **kw['plot'])
    elif command == 'plot-profile':
        for stream in iter_:
            channels = set(tr.stats.channel for tr in stream)
            for ch in channels:
                st2 = stream.select(channel=ch)
                fname = PLOT_PROFILE_FNAMES.format(root=path_out,
                                                   **st2[0].stats)
                _create_dir(fname)
                st2.plot_profile(fname, **kw['plot_profile'])
    elif command == 'stack':
        for stream in iter_:
            stack = stream.stack()
            write(stack, path_out, format, type='stack')
    elif command == 'profile':
        from rf.profile import get_profile_boxes, profile
        boxx = get_profile_boxes(**kw['boxes'])
        prof = profile(iter_, boxx, **kw['profile'])
        write(prof, path_out, format, type='profile')
    else:
        commands = [command] + list(commands)
        for stream in iter_:
            for command in commands:
                if command == 'data':
                    pass
                elif command == 'calc':
                    stream.rf(**kw['rf'])
                elif command == 'moveout':
                    stream.moveout(**kw['moveout'])
                else:
                    raise NotImplementedError
            write(stream, path_out, format)


def run_cli(args=None):
    """Command line interface of rf.

    After parsing call `run`.
    """
    from rf import __version__
    msg = ('*Note*: The command line tool is rather low priority in the '
           'development. It might also be depreciated in the future.')
    p = argparse.ArgumentParser(description=__doc__, epilog=msg)
    version = '%(prog)s ' + __version__
    p.add_argument('-v', '--version', action='version', version=version)
    msg = 'Configuration file to load (default: conf.json)'
    p.add_argument('-c', '--conf', default='conf.json', help=msg)

    sub = p.add_subparsers(title='commands', dest='command')
    msg = 'create config file in current directory'
    p_create = sub.add_parser('create', help=msg)
    msg = 'retrieve data for further processing'
    p_data = sub.add_parser('data', help=msg)
    msg = 'calculate receiver functions'
    p_calc = sub.add_parser('calc', help=msg)
    msg = 'perform move out correction'
    p_mout = sub.add_parser('moveout', help=msg)
    msg = 'stack receiver functions'
    p_stack = sub.add_parser('stack', help=msg)
    msg = 'stack receiver functions to profile'
    p_profile = sub.add_parser('profile', help=msg)
    msg = 'convert files to different format'
    p_conv = sub.add_parser('convert', help=msg)
    msg = 'print information about events, stations or waveform files'
    p_print = sub.add_parser('print', help=msg)
    msg = 'plot receiver functions'
    p_plot = sub.add_parser('plot', help=msg)
    msg = 'plot receiver function profile'
    p_plotp = sub.add_parser('plot-profile', help=msg)

    msg = 'create example files for tutorial'
    p_create.add_argument('-t', '--tutorial', help=msg, action='store_true')
    # the default='moveout' is an ugly work-around for
    # http://bugs.python.org/issue27227 and related issue9625
    msg = 'calculate receiver functions, perform moveout correction, optional'
    p_data.add_argument('commands', nargs='*', help=msg,
                        choices=('calc', 'moveout'), default='moveout')
    msg = 'perform also moveout correction'
    p_calc.add_argument('commands', nargs='*', help=msg,
                        choices=('moveout',), default='moveout')
    msg = "one of 'events', 'inventory' or filenames"
    p_print.add_argument('objects', nargs='+', help=msg)

    io = [p_calc, p_mout, p_conv, p_plot, p_stack, p_profile, p_plotp]
    for pp in io:
        msg = 'directory of files (SAC, Q) or basename of file (H5)'
        pp.add_argument('path_in', help=msg)
    io.append(p_data)
    for pp in io:
        msg = 'output directory or output file basename'
        pp.add_argument('path_out', help=msg)

    msg = 'new format (supported: Q, SAC or H5)'
    p_conv.add_argument('newformat', help=msg)

    msg = ('Use these flags to overwrite values in the config file. '
           'See the example configuration file for a description of '
           'these options.')
    g2 = p.add_argument_group('optional config arguments', description=msg)
    features_str = ('events', 'inventory', 'data', 'phase',
                    'moveout-phase', 'format')
    for f in features_str:
        g2.add_argument('--' + f, default=SUPPRESS)
    for f in DICT_OPTIONS:
        g2.add_argument('--' + f.replace('_', '-'), default=SUPPRESS)

    # Get command line arguments and start run
    args = vars(p.parse_args(args))
    # second part of the uggly work-around for argparse issue27227
    if args.get('commands') == 'moveout':
        args['commands'] = []
    try:
        run(**args)
    except ParseError as ex:
        p.error(ex)
