# -*- coding: utf-8 -*-
# Copyright 2013-2016 Tom Eulenfeld, MIT license
"""
Functions for receiver function profile calculation.
"""
import numpy as np
from rf.util import _add_processing_info, direct_geodetic


_LARGE_BOX_WIDTH = 2000


def _get_box(latlon0, azimuth, length, width=_LARGE_BOX_WIDTH, offset=0):
    """Create a single box."""
    from shapely.geometry import Polygon
    start = direct_geodetic(latlon0, azimuth, offset)
    azis = ((azimuth - 90) % 360, azimuth,
            (azimuth + 90) % 360, (azimuth + 180) % 360)
    dists = (width/2, length, width, length)
    latlon = start
    corners = []
    for a, d in zip(azis, dists):
        latlon = direct_geodetic(latlon, a, d)
        corners.append(latlon[::-1])
    box = {'poly': Polygon(corners),
           'length': length,
           'pos': offset + length/2,
           'latlon': direct_geodetic(start, azimuth, length/2)}
    return box


def get_profile_boxes(latlon0, azimuth, bins, width=_LARGE_BOX_WIDTH):
    """
    Create 2D boxes for usage in `profile()` function.

    :param tuple latlon0: coordinates of starting point of profile
    :param azimuth: azimuth of profile direction
    :param tuple bins: Edges of the distance bins in km (e.g. (0, 10, 20, 30))
    :param width: width of the boxes in km (default: large)
    :return: List of box dicts. Each box has the entries
        'poly' (shapely polygon with lonlat corners), 'length' (length in km),
        'pos' (midpoint of box in km from starting coordinates),
        'latlon' (midpoint of box as coordinates)
    """
    boxes = []
    for i in range(len(bins)-1):
        length = bins[i+1] - bins[i]
        box = _get_box(latlon0, azimuth, length, width, offset=bins[i])
        if i == 0:
            box['profile'] = {}
            box['profile']['latlon'] = latlon0
            box['profile']['azimuth'] = azimuth
            box['profile']['length'] = bins[-1] - bins[0]
            box['profile']['width'] = width
        boxes.append(box)
    return boxes


def _find_box(latlon, boxes, crs=None):
    """Return the box which encloses the coordinates."""
    import cartopy.crs as ccrs
    from shapely.geometry import Point
    if crs is None:
        latlons = [boxes[len(boxes)//2]['latlon']]
        latlon0 = np.median(latlons, axis=0)
        crs = ccrs.AzimuthalEquidistant(*latlon0[::-1])
    pc = ccrs.PlateCarree()
    p = crs.project_geometry(Point(*latlon[::-1]), pc)
    for box in boxes:
        poly = crs.project_geometry(box['poly'], pc)
        if p.within(poly):
            return box


@_add_processing_info
def profile(stream, boxes, crs=None):
    """
    Stack traces in stream by piercing point coordinates in defined boxes.

    :param stream: stream with pre-calculated piercing point coordinates
    :param boxes: boxes created with `get_profile_boxes()`
    :param crs: cartopy projection (default: AzimuthalEquidistant)
    :return: profile stream
    """
    stack = {}
    for tr in stream:
        ppoint = (tr.stats.pp_latitude, tr.stats.pp_longitude)
        box = _find_box(ppoint, boxes, crs=crs)
        if box is None:
            continue
        pos = box['pos']
        comp = tr.stats.channel[-1]
        key = (pos, comp)
        if key not in stack:
            header = {'box_pos': pos,
                      'box_length': box['length'],
                      'box_latitude': box['latlon'][0],
                      'box_longitude': box['latlon'][1],
                      'profile_latitude': boxes[0]['profile']['latlon'][0],
                      'profile_longitude': boxes[0]['profile']['latlon'][1],
                      'profile_azimuth': boxes[0]['profile']['azimuth'],
                      'profile_length': boxes[0]['profile']['length'],
                      'num': 1,
                      'sampling_rate': tr.stats.sampling_rate,
                      'channel': '??' + comp}
            for entry in ('slowness', 'phase', 'moveout', 'processing'):
                if entry in tr.stats:
                    header[entry] = tr.stats[entry]
            stack[key] = tr2 = tr.__class__(data=tr.data, header=header)
            if 'onset' in tr.stats:
                onset = tr.stats.onset - tr.stats.starttime
                tr2.stats.onset = tr2.stats.starttime + onset
        else:
            tr2 = stack[key]
            tr2.data = tr2.data + tr.data
            tr2.stats.num += 1
    for tr2 in stack.values():
        tr2.data = tr2.data / tr2.stats.num
    if hasattr(stream, 'iterable'):  # support tqdm objects
        cls = stream.iterable.__class__
    else:
        cls = stream.__class__
    try:
        profile = cls(traces=stack.values())
    except TypeError:  # stream can be an iterator
        from rf import RFStream
        profile = RFStream(traces=stack.values())
    profile.sort(['channel', 'box_pos'])
    profile.type = 'profile'
    return profile
