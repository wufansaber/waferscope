"""Complex whole-surface SPM plus a non-double-counted local FEM correction."""
import numpy as np
from numpy.polynomial.legendre import leggauss

from .far_field import radiation_amplitude
from .planar import background_fields
from .surface_scattering import SurfaceRadiation, radiation_kernel


def triangle_source(triangles, p, order=6):
    """Actual local triangular heights with the same incident/reflected reference field."""
    a, b, c = (triangles[:, i] for i in range(3))
    area2 = np.abs((b[:, 0]-a[:, 0])*(c[:, 1]-a[:, 1])-
                   (b[:, 1]-a[:, 1])*(c[:, 0]-a[:, 0]))
    nodes, weights = leggauss(order)
    nodes, weights = (nodes+1)/2, weights/2
    positions, values = [], []
    for u, wu in zip(nodes, weights):
        for v, wv in zip(nodes, weights):
            point = a+u*(b-a)+(1-u)*v*(c-a)
            positions.append(point[:, :2])
            values.append(point[:, 2]*area2*(1-u)*wu*wv)
    positions, values = np.concatenate(positions), np.concatenate(values)
    plane = np.column_stack((positions, np.zeros(len(positions))))
    fields, _ = background_fields(plane, p)
    return positions, values[:, None]*fields


def triangle_amplitude(positions, weights, directions, p):
    directions = np.asarray(directions)
    shape = directions.shape[:-1]
    flat = directions.reshape((-1, 3))
    result = np.empty((len(flat), 2), complex)
    for start in range(0, len(flat), 64):
        q = flat[start:start+64]
        transform = np.exp(-1j*p.k0*(q[:, :2]@positions.T))@weights
        result[start:start+64] = np.einsum('nij,nj->ni', radiation_kernel(q, p), transform)
    return result.reshape(shape+(2,))


def combine(whole_spm, local_spm, local_fem_with, local_fem_without):
    common = whole_spm-local_spm
    return common+local_fem_with, common+local_fem_without


class HybridRadiation:
    def __init__(self, data, p):
        self.data, self.p = data, p
        axis = data['whole_axis_um']
        self.whole = SurfaceRadiation(axis, axis, data['whole_weighted_fields'], p, padding=4)
        self._directions = self._pair = None

    def components(self, directions):
        directions = np.asarray(directions)
        flat = directions.reshape((-1, 3))
        shape = directions.shape[:-1]+(2,)
        whole = self.whole.amplitude(directions)
        local = triangle_amplitude(self.data['local_spm_positions_um'],
                                   self.data['local_spm_weights'], directions, self.p)
        fields = [radiation_amplitude(self.data['source_points_'+name],
            self.data['source_polarization_weights_'+name], flat, self.p).reshape(shape)
            for name in ('with', 'without')]
        with_particle, without = combine(whole, local, *fields)
        return dict(with_particle=with_particle, without=without, whole_spm=whole,
                    local_spm=local, local_fem_with=fields[0], local_fem_without=fields[1])

    def amplitude(self, directions, state='with'):
        if state not in ('with', 'without'):
            raise ValueError('Hybrid state must be with or without particle')
        if self._directions is None or not np.array_equal(directions, self._directions):
            values = self.components(directions)
            self._directions = np.array(directions, copy=True)
            self._pair = (values['with_particle'], values['without'])
        return self._pair[0 if state == 'with' else 1]
