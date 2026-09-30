"""Native orthogonal near-field samples, with a legacy YZ archive adapter."""
from dataclasses import dataclass

import numpy as np


PLANE_AXES = {'xz': (0, 2, 1), 'xy': (0, 1, 2), 'yz': (1, 2, 0)}


@dataclass(frozen=True)
class NearPlane:
    name: str
    u: np.ndarray
    v: np.ndarray
    offset: float

    @property
    def points(self):
        horizontal, vertical, fixed = PLANE_AXES[self.name]
        uu, vv = np.meshgrid(self.u, self.v)
        points = np.empty((uu.size, 3))
        points[:, horizontal] = uu.ravel()
        points[:, vertical] = vv.ravel()
        points[:, fixed] = self.offset
        return points

    def metadata(self):
        horizontal, vertical, fixed = PLANE_AXES[self.name]
        return dict(plane=self.name.upper(), horizontal_axis='XYZ'[horizontal],
                    vertical_axis='XYZ'[vertical], fixed_axis='XYZ'[fixed],
                    position_um=float(self.offset), shape=[len(self.v), len(self.u)])


def make_planes(p, samples=101):
    p.validate()
    if not isinstance(samples, int) or samples < 3:
        raise ValueError('Near-field samples must be an integer >= 3')
    radius = p.domain_radius_um
    lateral = np.linspace(-.6*radius, .6*radius, samples)
    height = np.linspace(-.299*radius, .75*radius, samples)
    center = (p.diameter_nm/2+p.gap_nm)/1000
    # Keep even the XY rectangle corners inside the physical sphere, outside PML.
    half_width = min(.6*radius, .98*np.sqrt((radius**2-center**2)/2))
    xy_axis = np.linspace(-half_width, half_width, samples)
    return (NearPlane('xz', lateral, height, 0.),
            NearPlane('xy', xy_axis, xy_axis, center),
            NearPlane('yz', lateral, height, 0.))


def sample_planes(mesh, field, p, planes, state, data):
    from .solver import evaluate_fields
    for plane in planes:
        prefix = 'near_'+plane.name+'_'
        data[prefix+'u_um'], data[prefix+'v_um'] = plane.u, plane.v
        data[prefix+'position_um'] = np.asarray(plane.offset)
        electric, magnetic, _ = evaluate_fields(mesh, field, plane.points, p)
        shape = (len(plane.v), len(plane.u), 3)
        data[prefix+'e_'+state] = electric.reshape(shape)
        data[prefix+'hz0_'+state] = magnetic.reshape(shape)


def finish_planes(data, planes):
    for plane in planes:
        for quantity in ('e', 'hz0'):
            prefix = 'near_'+plane.name+'_'+quantity+'_'
            data[prefix+'change'] = data[prefix+'with']-data[prefix+'without']
    # Existing audits and saved-result readers still use these YZ-only names.
    data['near_y_um'], data['near_z_um'] = data['near_yz_u_um'], data['near_yz_v_um']
    for quantity in ('e', 'hz0'):
        for state in ('with', 'without', 'change'):
            data['near_'+quantity+'_'+state] = data['near_yz_'+quantity+'_'+state]


def saved_planes(data):
    planes = []
    for name in PLANE_AXES:
        prefix = 'near_'+name+'_'
        if prefix+'e_with' in data:
            planes.append(NearPlane(name, data[prefix+'u_um'], data[prefix+'v_um'],
                                    float(data[prefix+'position_um'])))
    if not planes:
        planes.append(NearPlane('yz', data['near_y_um'], data['near_z_um'], 0.))
    return planes


def electric_field(data, plane):
    key = 'near_'+plane.name+'_e_with'
    return data[key] if key in data else data['near_e_with']
