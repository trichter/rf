# -*- coding: utf-8 -*-
"""
Functions for receiver function plotting.
"""

import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator, MaxNLocator
import numpy as np
import warnings


def plot_rf(stream, fname=None, norm=1., fig_width=7., trace_height=0.5,
            stack_height=0.5,
            fill=False, window=None, downsample=None, title=True,
            info=[('back_azimuth', u'baz (°)', 'b'),
                  ('distance', u'dist (°)', 'r')]):
    """
    Create receiver function plot.

    :param fname: Filename to save plot to. Can be None. In this case
        the figure is left open.
    :param fig_width: Width of figure in inches.
    :param trace_height: Height of one trace in inches.
    :param fill: Waether to fill receiver functions or not.
    :param downsample: Downsample to frequency (in Hz) with
        Stream.decimate. Filtering is not performed. When saving in a
        vector format the plot size can be reduced in this way.
    :param title: Print seed id as a title.
    :param info: Plot one additional axes showing maximal two entries of
        the stats object. Each entry in this list is a list consisting of
        three entries: key, label and color.
        info can be None. In this case no additional axes is plotted.
    """

    if len(stream) == 0:
        return
    if window:
        for tr in stream:
            tr.trim(tr.stats.onset + window[0], tr.stats.onset + window[1])
    if downsample:
        for tr in stream:
            tr.decimate(int(round(tr.stats.sampling_rate)) // downsample,
                        no_filter=True)
    # calculate lag times
    stats = stream[0].stats
    N = len(stream)
    t0 = stats.onset - stats.starttime
    t2 = stats.endtime - stats.starttime
    times = np.linspace(-t0, t2 - t0, stats.npts, endpoint=True)
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
    ax2 = fig.add_axes([fl, 1 - ft - hs, fw2, hs], sharex=ax1)
    if info:
        ax3 = fig.add_axes(
            [1 - fr - fw3, fb, fw3, h * (N + 2)], sharey=ax1)
        info = list(info)
        info[0] = [ax3] + list(info[0])
        if len(info) > 1:
            ax4 = ax3.twiny()
            info[1] = [ax4] + list(info[1])
    # plot stack and individual receiver functions
    stack = stream.stack()
    if len(stack) > 1:
        warnings.warn('Different stations in one RF plot.')

    def _rf_fill(ax, t, d, i):
        ax.fill_between(t, d + i, i, where=d >= 0, lw=0., facecolor='k')
        ax.fill_between(t, d + i, i, where=d < 0, lw=0, facecolor='grey')

    def _plot(ax, t, d, i):
        if fill:
            _rf_fill(ax, t, d, i)
            ax.plot(t, d + i, 'k')
        else:
            ax.plot(t, d + i, 'k')
    _plot(ax2, times, stack[0].data, 0)
    for i, tr in enumerate(stream):
        _plot(ax1, times, tr.data * norm, i + 1)
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
    for l in ax2.get_xticklabels():
        l.set_visible(False)
    ax2.yaxis.set_major_locator(MaxNLocator(4))
    for l in ax2.get_yticklabels():
        l.set_fontsize('small')
    # plot title and save plot
    if title:
        bbox = dict(boxstyle='round', facecolor='white', alpha=0.8, lw=0)
        text = '%s traces  %s' % (len(stream), stack[0].id)
        ax2.annotate(text, (1 - 0.5 * fr, 1 - 0.5 * ft),
                     xycoords='figure fraction', va='top', ha='right',
                     bbox=bbox, clip_on=False)
    if fname:
        fig.savefig(fname)
        plt.close(fig)
