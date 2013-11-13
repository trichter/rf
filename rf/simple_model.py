"""
Simple move out and piercing point calculation.

This module will probably be replaced by something more sophisticated in
future.
"""
import os.path
from math import floor
import numpy as np
from obspy.core.util.geodetics import kilometer2degrees
try:
    from geographiclib.geodesic import Geodesic
except ImportError:
    import warnings
    msg = 'Geographiclib import error. Ppoint calculation will not work.'
    warnings.warn(msg)


_model_cache = {}


def load_model(fname='iasp91'):
    """
    Load model from file.

    :param fname: path to model file or 'iasp91'
    :return: :class:`~rf.simple.model.SimpleModel` instance
    """
    try:
        return _model_cache[fname]
    except KeyError:
        pass
    fname_key = fname
    if fname == 'iasp91':
        fname = os.path.join(os.path.dirname(__file__), 'data', 'iasp91.dat')
    values = np.loadtxt(fname, unpack=True)
    try:
        z, vp, vs, n = values
        n.astype(int)
    except ValueError:
        n = None
        z, vp, vs = values
    _model_cache[fname_key] = model = SimpleModel(z, vp, vs, n)
    return model


def _interpolate_n(val, n):
    vals = [np.linspace(val[i], val[i + 1], n[i + 1] + 1, endpoint=False)
            for i in range(len(val) - 1)]
    return np.hstack(vals + [np.array([val[-1]])])


class SimpleModel(object):
    """
    Simple 1D velocity model.
    """
    def __init__(self, z, vp, vs, n=None):
        assert len(z) == len(vp) == len(vs)
        if n is not None:
            z = _interpolate_n(z, n)
            vp = _interpolate_n(vp, n)
            vs = _interpolate_n(vs, n)
        self.z = z[:-1]
        self.dz = np.diff(z)
        self.vp = vp[:-1]
        self.vs = vs[:-1]
        self.t_ref = {}

    def _calc_vertical_slowness(self, slowness, phase='PS'):
        """
        Calculate vertical slowness of P and S wave.
        """
        qp, qs = 0, 0
        if 'P' in phase:
            qp = np.sqrt(self.vp ** (-2) - slowness ** 2)
        if 'S' in phase:
            qs = np.sqrt(self.vs ** (-2) - slowness ** 2)
        return qp, qs

    def _calc_phase_delay(self, slowness, phase='PS'):
        """
        Calculate travel time delay between direct wave and converted phase.
        """
        qp, qs = self._calc_vertical_slowness(slowness, phase=phase)
        dt = (qp * phase.count('P') + qs * phase.count('S') -
              2 * (qp if phase[0] == 'P' else qs)) * self.dz
        return np.cumsum(dt)

    def calculate_phase_delays(self, slowness, phase='Ps', ref=6.4):
        """
        Calculate travel time delays for slowness and reference slowness.

        Travel time delays are calculated between the direct wave and the
        converted phase or multiples.

        :param slowness: ray parameter in s/deg
        :param phase: 'Ps', 'Sp' or multiples
        :param ref: reference ray parameter in s/deg
        :return: original times, times stretched to reference slowness
        """
        assert len(phase) % 2 == 0
        slowness = slowness * kilometer2degrees(1)
        ref = ref * kilometer2degrees(1)
        phase = phase.upper()
        try:
            t_ref = self.t_ref[phase]
        except KeyError:
            self.t_ref[phase] = t_ref = self._calc_phase_delay(ref, phase)
        t = self._calc_phase_delay(slowness, phase)
        if phase[0] == 'S':
            t_ref = -t_ref
            t = -t
        index = np.nonzero(np.isnan(t))[0][0] - 1
        t = t[:index]
        t_ref = t_ref[:index]
        return (np.hstack((-t[1:10][::-1], t)),
                np.hstack((-t_ref[1:10][::-1], t_ref)))

    def moveout(self, stream, phase='Ps', ref=6.4):
        """
        In-place moveout correction.

        :param stream: stream with stats attributes onset and slowness.
        :param phase: 'Ps', 'Sp', 'Ppss' or other multiples
        :param ref: reference ray parameter in s/deg
        """
        for tr in stream:
            st = tr.stats
            index0 = int(floor((st.onset - st.starttime) * st.sampling_rate))
            t0, t1 = self.calculate_phase_delays(st.slowness, phase=phase,
                                                 ref=ref)
            S_multiple = phase[0].upper() == 'S' and len(phase) > 3
            if S_multiple:
                time0 = st.starttime - st.onset + index0 * st.delta
                old_data = tr.data[:index0][::-1]
                t = -time0 - np.arange(index0) * st.delta
                new_t = -np.interp(-t, -t0, -t1, left=0, right=0)
                data = np.interp(-t, -new_t, old_data, left=0., right=0.)
                tr.data[:index0] = data[::-1]
            else:
                if t0[-1] > t1[-1]:
                    index0 += 1
                time0 = st.starttime - st.onset + index0 * st.delta
                old_data = tr.data[index0:]
                t = time0 + np.arange(len(tr) - index0) * st.delta
                # stretch old times to new times
                new_t = np.interp(t, t0, t1, left=0, right=0)
                # interpolate data at new times to data samples
                data = np.interp(t, new_t, old_data, left=0., right=0.)
                tr.data[index0:] = data

    def ppoint_distance(self, depth, slowness, phase='S'):
        """
        Calculate horizontal distance of piercing point to station.

        :param depth: depth of interface in km
        :param slowness: ray parameter in s/deg
        :param phase: 'P' or 'S' for P wave or S wave. Multiples possible.
        :return: horizontal distance in km
        """
        assert len(phase) % 2 == 1
        phase = phase.upper()
        slowness = slowness * kilometer2degrees(1)
        xp, xs = 0., 0.
        qp, qs = self._calc_vertical_slowness(slowness, phase=phase)
        if 'P' in phase:
            xp = np.cumsum(self.dz * slowness / qp)
        if 'S' in phase:
            xs = np.cumsum(self.dz * slowness / qs)
        x = xp * phase.count('P') + xs * phase.count('S')
        z = self.z
        index = np.nonzero(depth < z)[0][0] - 1
        return x[index] + ((x[index + 1] - x[index]) *
                           (depth - z[index]) / (z[index + 1] - z[index]))

    def ppoint(self, stats, depth, phase='S'):
        """
        Piercing point calculation.

        Piercing point coordinates are saved in the plat and plon attributes of
        the stats objects.

        :param stats: Stats object with attributes slowness, back_azimuth and
            station coordinates.
        :param depth: depth of interface in km
        :param phase: 'P' for piercing points of P wave, 'S' for piercing
            points of S wave. Multiples are possible, too.
        :return: piercing point latitude and longitude
        """
        dr = self.ppoint_distance(depth, stats.slowness, phase=phase)
        lat = stats.station_latitude
        lon = stats.station_longitude
        az = stats.back_azimuth
        result = Geodesic.WGS84.Direct(lat, lon, az, 1000 * dr)
        stats.plat = result['lat2']
        stats.plon = result['lon2']
        return stats.plat, stats.plon
