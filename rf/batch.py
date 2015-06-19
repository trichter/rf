"""
rf batch command line utility
"""

import argparse
import os
from os.path import join
import shutil
import sys
from obspy import read_inventory, readEvents
from rf.rfstream import read_rf, rfstats, RFStream, set_index
from pkg_resources import resource_filename
try:
    from progressbar import ProgressBar
except ImportError:
    ProgressBar = None

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
    inventory = read_inventory(inventory)
    channels = inventory.get_contents()['channels']
    stations = list(set(ch.rsplit('.', 1)[0] for ch in channels))
    one_channel = {ch.rsplit('.', 1)[0]: ch for ch in channels}
    if events is not None:
        events = readEvents(events)
        yield len(stations) * len(events)
        for event in events:
            for station in stations:
                seed_id = one_channel[station][:-1] + '?'
                net, sta, loc, cha = seed_id.split('.')
                stats = {'network': net, 'station': sta, 'location': loc,
                         'channel': cha}
                if rf:
                    stats['event'] = event
                    stats['seed_id'] = seed_id
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


def batch_rf(init, events, inventory, path, format='H5',
             request_window=None, phase='P', dist_range=None, **rf_kwargs):
    get_waveform = init()
    root = phase + 'rf'
    _check_path(join(path, root))
    method = phase[-1]
    if dist_range is None:
        dist_range = (30, 90) if method == 'P' else (60, 85)
    if request_window is None:
        request_window = (-50, 150) if method == 'P' else (-80, 50)
    for kwargs, event, coords in _iter(events, inventory, rf=True):
        stats = rfstats(station=coords, event=event,
                        phase=phase, dist_range=dist_range)
        if not stats:
            continue
        kwargs.update({'starttime': stats.onset + request_window[0],
                       'endtime': stats.onset + request_window[1]})
        stream = get_waveform(**kwargs)
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


def batch_moveout(phase, events, inventory, path, root, format,
                  kwargs_moveout, **kwargs):
    if root is None:
        root = phase[0] + 'rf'
    for stats in _iter(events, inventory):
        stream = _read(stats, path, root, format)
        if stream is None:
            continue
        stream.moveout(phase, **kwargs_moveout)
        if len(stream) > 0:
            _write(stream, path, '%s_%s' % (root, phase), format)


def batch_convert(events, inventory, path, root, newformat, format, **kwargs):
    for stats in _iter(events, inventory):
        stream = _read(stats, path, root, format)
        if stream is None:
            continue
        _write(stream, path, join(newformat, root), newformat)


def batch_plot(events, inventory, path, root, format, kwargs_plot, **kwargs):
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


def batch_stack(inventory, path, root, format, **kwargs):
    for stats in _iter(None, inventory):
        stream = _read(stats, path, root, format)
        stream.stack()
        _write(stream, path, 'stack_%s' % root, format, stack=True)

CONFIG = ['events', 'inventory', 'init', 'request_window',
          'dist_range', 'path', 'format',
          'phase', 'filter', 'window', 'downsample', 'rotate',
          'deconvolve', 'source_component', 'winsrc']
CONFIG_FREQ = ['water', 'gauss', 'tschift']
CONFIG_TIME = ['winrsp', 'winrf', 'spiking']

def _read_config(fname, conf={}):
    with open(fname) as f:
        exec(compile(f.read(), fname, 'exec'), conf)
    return conf

def _slice_config(conf):
    if conf['deconvolve'] == 'freq':
        CONFIG.extend(CONFIG_FREQ)
    else:
        CONFIG.extend(CONFIG_TIME)
    return {k: conf[k] for k in CONFIG if k in conf}


def main(args=None):
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(title='subcommands', dest='subcommand')
    parser_init = subparsers.add_parser('init', help='initialize new project')
    parser_init.add_argument('path', help='path for new rf project')

    help = 'calculate receiver functions'
    parser_calc = subparsers.add_parser('calc', help=help)
    help = 'check events, station and configuration file'
    parser_check = subparsers.add_parser('check', help=help)
    help = 'move out correction'
    parser_mout = subparsers.add_parser('moveout', help=help)
    help = 'convert files to different format'
    parser_conv = subparsers.add_parser('convert', help=help)
    help = 'plot receiver functions'
    parser_plot = subparsers.add_parser('plot', help=help)
    help = 'stack receiver functions'
    parser_stack = subparsers.add_parser('stack', help=help)
    help = 'init project with example files for a tutorial'
    parser_init.add_argument('-t', '--tutorial', help=help,
                             action='store_true')
    help = 'path with config file name (default: config.py)'
    for p in [parser_calc, parser_check, parser_mout, parser_conv, parser_plot,
              parser_stack]:
        p.add_argument('-c', '--conf', default='conf.py', help=help)
    for p in [parser_mout, parser_conv, parser_plot, parser_stack]:
        p.add_argument('root', help='directory or basename of file')
    help = ('Phase (P, S, PP, etc.). The last letter (P or S) is used as a'
            'specification of the used method (P or S receiver function).'
            'E.g. PP means use PP phase for P-receiver function calculation')
    parser_calc.add_argument('phase', help=help)
    parser_mout.add_argument('phase', help='phase to align (Ps, Ppps, etc.)')
    help = 'convert to other format (supported: Q, SAC or H5)'
    parser_conv.add_argument('format', help=help)

    args = parser.parse_args(args)
    if args.subcommand == 'init':
        path = args.path
        _check_path(path)
        os.mkdir(path)
        fnames = ['conf.py']
        if args.tutorial:
            fnames.extend(['example_events.xml', 'example_inventory.xml',
                           'example_data.mseed'])
        for fname in fnames:
            src = resource_filename('rf', 'example/%s' % fname)
            dst = os.path.join(path, fname)
            shutil.copyfile(src, dst)
        print('New rf project initialized in directory %s' % path)
    else:
        if args.subcommand == 'calc':
            conf = {'phase': args.phase, 'method': args.phase[-1].upper()}
        else:
            conf = {'phase': 'P', 'method': 'P'}
        conf = _read_config(args.conf, conf)
        if args.subcommand == 'calc':
            conf = _slice_config(conf)
            batch_rf(**conf)
        elif args.subcommand == 'check':
            print('Check %s' % readEvents(conf['events']).__str__(True))
            print('Check %s' % read_inventory(conf['inventory']))
        elif args.subcommand == 'moveout':
            conf.pop('phase')
            batch_moveout(phase=args.phase, root=args.root, **conf)
        elif args.subcommand == 'convert':
            batch_convert(root=args.root, newformat=args.format, **conf)
        elif args.subcommand == 'plot':
            batch_plot(root=args.root, **conf)
        elif args.subcommand == 'stack':
            batch_stack(root=args.root, **conf)
        else:
            print('Subcommand not availlable')


if __name__ == '__main__':
    main()
