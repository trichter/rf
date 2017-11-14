# -*- coding: utf-8 -*-
# Copyright 2013-2016 Tom Eulenfeld, MIT license
"""
Functions for receiver function plotting.
"""
import warnings

import matplotlib.patheffects as PathEffects
from matplotlib.ticker import (AutoMinorLocator, FixedLocator, FixedFormatter,
                               MaxNLocator)
import matplotlib.pyplot as plt
import numpy as np


def plot_rf(stream, fname=None, fig_width=7., trace_height=0.5,
            stack_height=0.5, scale=1, fillcolors=(None, None), trim=None,
            info=(('back_azimuth', u'baz (°)', 'C0'),
                  ('distance', u'dist (°)', 'C3'))):
    """
    Plot receiver functions.

    :param stream: stream to plot
    :param fname: filename to save plot to. Can be None. In this case
        the figure is left open.
    :param fig_width: width of figure in inches
    :param trace_height: height of one trace in inches
    :param stack_height: height of stack axes in inches
    :param scale: scale for individual traces
    :param fillcolors: fill colors for positive and negative wiggles
    :param trim: trim stream relative to onset before plotting using
         `~.rfstream.RFStream.slice2()`
    :param info: Plot one additional axes showing maximal two entries of
        the stats object. Each entry in this list is a list consisting of
        three entries: key, label and color.
        info can be None. In this case no additional axes is plotted.
    """

    if len(stream) == 0:
        return
    if trim:
        stream = stream.slice2(*trim, reftime='onset')
    N = len(stream)
    # calculate lag times
    stats = stream[0].stats
    times = stream[0].times() - (stats.onset - stats.starttime)
    # calculate axes and figure dimensions
    # big letters: inches, small letters: figure fraction
    H = trace_height
    HS = stack_height
    FB = 0.5
    FT = 0.2
    DW = 0.1
    FH = H * (N + 2) + HS + FB + FT + DW
    h = H / FH
    hs = HS / FH
    fb = FB / FH
    ft = FT / FH
    FL = 0.5
    FR = 0.2
    FW = fig_width
    FW3 = 0.8
    FW2 = FW - FL - FR - (DW + FW3) * bool(info)
    fl = FL / FW
    fr = FR / FW
    fw2 = FW2 / FW
    fw3 = FW3 / FW
    # init figure and axes
    fig = plt.figure(figsize=(FW, FH))
    ax1 = fig.add_axes([fl, fb, fw2, h * (N + 2)])
    if info:
        ax3 = fig.add_axes(
            [1 - fr - fw3, fb, fw3, h * (N + 2)], sharey=ax1)
        info = list(info)
        info[0] = [ax3] + list(info[0])
        if len(info) > 1:
            ax4 = ax3.twiny()
            info[1] = [ax4] + list(info[1])
    # plot individual receiver functions

    def _plot(ax, t, d, i):
        c1, c2 = fillcolors
        if c1:
            ax.fill_between(t, d + i, i, where=d >= 0, lw=0., facecolor=c1)
        if c2:
            ax.fill_between(t, d + i, i, where=d < 0, lw=0., facecolor=c2)
        ax.plot(t, d + i, 'k')
    max_ = max(np.max(np.abs(tr.data)) for tr in stream)
    for i, tr in enumerate(stream):
        _plot(ax1, times, tr.data / max_ * scale, i + 1)
    # plot right axes with header information
    for ax, header, label, color in info:
        data = [tr.stats[header] for tr in stream]
        ax.plot(data, 1 + np.arange(len(stream)), '.' + color, mec=color)
        ax.set_xlabel(label, color=color, size='small')
        if header == 'back_azimuth':
            ax.set_xticks(np.arange(5) * 90)
            ax.set_xticklabels(['0', '', '180', '', '360'], size='small')
        else:
            ax.xaxis.set_major_locator(MaxNLocator(4))
            for l in ax.get_xticklabels():
                l.set_fontsize('small')
        ax.xaxis.set_minor_locator(AutoMinorLocator())
    # set x and y limits
    ax1.set_xlim(times[0], times[-1])
    ax1.set_ylim(-0.5, N + 1.5)
    ax1.set_yticklabels('')
    ax1.set_xlabel('time (s)')
    ax1.xaxis.set_minor_locator(AutoMinorLocator())

    # plot stack
    stack = stream.stack()
    if len(stack) > 1:
        warnings.warn('Different stations or channels in one RF plot.')
    elif len(stack) == 1:
        ax2 = fig.add_axes([fl, 1 - ft - hs, fw2, hs], sharex=ax1)
        _plot(ax2, times, stack[0].data, 0)
        for l in ax2.get_xticklabels():
            l.set_visible(False)
        ax2.yaxis.set_major_locator(MaxNLocator(4))
        for l in ax2.get_yticklabels():
            l.set_fontsize('small')
        # annotate plot with seed id
        bbox = dict(boxstyle='round', facecolor='white', alpha=0.8, lw=0)
        text = '%s traces  %s' % (len(stream), stack[0].id)
        ax2.annotate(text, (1 - 0.5 * fr, 1 - 0.5 * ft),
                     xycoords='figure fraction', va='top', ha='right',
                     bbox=bbox, clip_on=False)
    # save plot
    if fname:
        fig.savefig(fname)
        plt.close(fig)
    else:
        return fig


def _get_geoaxes(crs=None, latlons=None):
    """Return cartopy geoaxis"""
    if crs is None:
        from cartopy.crs import AzimuthalEquidistant
        latlon0 = np.median(latlons, axis=0)
        crs = AzimuthalEquidistant(*latlon0[::-1])
    return plt.axes(projection=crs)


def __pc():
    from cartopy.crs import PlateCarree as PC
    return PC()


def plot_stations(inventory, label_stations=True, ax=None, crs=None, **kwargs):
    """
    Plot stations.

    :param inventory: station inventory
    :param label_stations: weather to label stations
    :param ax: geoaxes (default None: new ax will be created)
    :param crs: coordinate reference system for new geoaxis, (default: None,
        then AzimuthalEquidistant projection with appropriate center is used.)
    :param \*\*kwargs: other kwargs are passed to ax.scatter() call
    """
    try:
        # assume inventory to be stream
        data = [((tr.stats.station_latitude, tr.stats.station_longitude),
                 tr.stats.station) for tr in inventory]
        latlons, names = zip(*list(set(data)))
    except AttributeError:
        # inventory is Inventory
        latlons, names = zip(*[((sta.latitude, sta.longitude), sta.code)
                               for net in inventory for sta in net])
    if ax is None:
        ax = _get_geoaxes(crs=crs, latlons=latlons)
    kw = dict(s=200, marker='v', c='darkred', linewidth=0.5, zorder=3)
    kw.update(kwargs)
    ax.scatter(*list(zip(*latlons))[::-1], transform=__pc(), **kw)
    if label_stations:
        path_effect = PathEffects.withStroke(linewidth=3, foreground="white")
        kw = {'xycoords': __pc()._as_mpl_transform(ax),
              'xytext': (10, 0), 'textcoords': 'offset points', 'zorder': 4,
              'path_effects': [path_effect]}
        for latlon, name in zip(latlons, names):
            ax.annotate(name, latlon[::-1], **kw)
    return ax


def plot_ppoints(ppoints, inventory=None, label_stations=True, ax=None,
                 crs=None, **kwargs):
    """
    Plot piercing points with stations.

    :param ppoints: list of (lat, lon) tuples of piercing points
    :param inventory, label_stations: plot stations, see `plot_stations`
    :param ax: geoaxes (default None: new ax will be created)
    :param crs: coordinate reference system for new geoaxis, (default: None,
        then AzimuthalEquidistant projection with appropriate center is used.)
    :param \*\*kwargs: other kwargs are passed to ax.scatter() call
    """
    if ax is None:
        ax = _get_geoaxes(crs=crs, latlons=ppoints)
    if inventory is not None:
        plot_stations(inventory, label_stations=label_stations, ax=ax)
    kw = dict(s=50, marker='x', color='k', alpha=0.2, zorder=2)
    kw.update(kwargs)
    ax.scatter(*list(zip(*ppoints))[::-1], transform=__pc(), **kw)
    return ax


def plot_profile_map(boxes, inventory=None, label_stations=True, ppoints=None,
                     ax=None, crs=None, **kwargs):
    """
    Plot profile map with stations and piercing points.

    :param boxes: boxes created with `~.profile.get_profile_boxes()`
    :param inventory, label_stations: plot stations, see `plot_stations`
    :param ppoints: list of (lat, lon) tuples of piercing points,
        see `plot_ppoints`
    :param ax: geoaxes (default None: new ax will be created)
    :param crs: coordinate reference system for new geoaxis, (default: None,
        then AzimuthalEquidistant projection with appropriate center is used.)
    :param \*\*kwargs: other kwargs are passed to ax.add_geometries() call
    """
    if ax is None:
        latlons = [boxes[len(boxes)//2]['latlon']]
        ax = _get_geoaxes(crs=crs, latlons=latlons)
    if inventory is not None:
        plot_stations(inventory, label_stations=label_stations, ax=ax)
    if ppoints is not None:
        plot_ppoints(ppoints, ax=ax)
    kw = dict(facecolor='none', edgecolor='0.8', zorder=1)
    kw.update(kwargs)
    for box in boxes:
        ax.add_geometries([box['poly']], crs=__pc(), **kw)
    return ax


def plot_profile(profile, fname=None, scale=1, fillcolors=('C3', 'C0'),
                 trim=None, top=None, moveout_model='iasp91'):
    """
    Plot receiver function profile.

    :param profile: stream holding the profile
    :param fname: filename to save plot to. Can be None. In this case
        the figure is left open.
    :param scale: scale for individual traces
    :param fillcolors: fill colors for positive and negative wiggles
    :param trim: trim stream relative to onset before plotting using
         `~.rfstream.RFStream.slice2()`
    :param top: show second axes on top of profile with additional information.
        Valid values: 'hist' - Plot histogram showing the number of receiver
        functions stacked in the corresponding bin
    :param moveout_model: string with model filename. Will be loaded into a
        `~.simple_model.SimpleModel` object to calculate depths for
        tick labels.
    """
    if len(profile) == 0:
        return
    if trim:
        profile = profile.slice2(*trim, reftime='onset')
    fig = plt.figure()
    ax = fig.add_axes([0.1, 0.1, 0.8, 0.7])
    widths = [tr.stats.box_length for tr in profile]
    pad = max(1, scale) * min(widths)
    xlim = (min(tr.stats.box_pos for tr in profile) - pad,
            max(tr.stats.box_pos for tr in profile) + pad)
    max_ = max(np.max(np.abs(tr.data)) for tr in profile)
    for tr in profile:
        x = tr.stats.box_pos + scale * tr.data / max_ * min(widths)
        y = tr.times() - (tr.stats.onset - tr.stats.starttime)
        ax.plot(x, y, 'k')
        c1, c2 = fillcolors
        if c1:
            ax.fill_betweenx(y, x, tr.stats.box_pos,
                             where=x >= tr.stats.box_pos, facecolor=c1)
        if c2:
            ax.fill_betweenx(y, x, tr.stats.box_pos,
                             where=x < tr.stats.box_pos, facecolor=c2)
    ax.set_xlabel('distance (km)')
    ax.set_ylim(max(y), min(y))
    ax.set_ylabel('time (s)')
    if moveout_model:
        from rf.simple_model import load_model
        model = load_model(moveout_model)
        phase = profile[0].stats.moveout
        slowness = profile[0].stats.slowness
        pd = model.calculate_delay_times(phase=phase, slowness=slowness)
        ax2 = ax.twinx()
        ax.get_shared_y_axes().join(ax, ax2)
        dkm = 50
        if profile[0].stats.endtime - profile[0].stats.onset > 50:
            dkm = 200
        d1 = np.arange(20) * dkm
        d2 = np.arange(100) * dkm / 5
        t1 = np.interp(d1, model.z, pd)
        t2 = np.interp(d2, model.z, pd)
        myLocator = FixedLocator(t1)
        myMinorLocator = FixedLocator(t2)
        myFormatter = FixedFormatter([str(i) for i in d1])
        ax2.yaxis.set_major_locator(myLocator)
        ax2.yaxis.set_minor_locator(myMinorLocator)
        ax2.yaxis.set_major_formatter(myFormatter)
        ax2.set_ylabel('depth (km)')
        ax2.set_ylim(ax.get_ylim())
    if top is not None:
        ax3 = fig.add_axes([0.1, 0.85, 0.8, 0.1], sharex=ax)
    if top == 'hist':
        left = [tr.stats.box_pos - tr.stats.box_length / 2 for tr in profile]
        height = [tr.stats.num for tr in profile]
        ax3.bar(left, height, widths, color='cadetblue')
        plt.setp(ax3.get_xticklabels(), visible=False)
        ax3.spines['top'].set_color('none')
        ax3.spines['right'].set_color('none')
        ax3.spines['left'].set_color('none')
        ax3.xaxis.set_ticks_position('bottom')
        ax3.yaxis.set_ticks_position('left')
        ax3.set_yticks(ax3.get_ylim())
    elif top is not None:
        raise NotImplementedError("'%s' not supported for top parameter" % top)
    ax.set_xlim(*xlim)
    if fname:
        fig.savefig(fname)
        plt.close(fig)
    else:
        return fig
