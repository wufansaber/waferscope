"""Finite trench motif surrounded by a flat wafer within an S4 supercell."""
from dataclasses import dataclass
import numpy as np


@dataclass(frozen=True)
class FinitePatch:
    length_um: float = .8
    trench_count: int = 2

    def centers(self, scene):
        return (np.arange(self.trench_count)-(self.trench_count-1)/2)*scene.pitch_um

    def validate(self, scene):
        if (not np.isfinite(self.length_um) or self.length_um <= 0
                or isinstance(self.trench_count, bool) or not isinstance(self.trench_count, int)
                or not 2 <= self.trench_count <= 32 or self.trench_count % 2):
            raise ValueError('Positive finite length and an even number of 2..32 trenches required')
        width = scene.pitch_um-scene.ridge_width_um
        if scene.kind != 'grating' or width <= 0:
            raise ValueError('Finite etched trenches require positive trench width')
        extent = max(abs(self.centers(scene)))+width/2
        if self.length_um >= scene.lx_um or 2*extent >= scene.ly_um:
            raise ValueError('Finite trenches must have flat buffer regions on all lateral sides')
        if scene.diameter_um >= min(scene.ridge_width_um, self.length_um):
            raise ValueError('The centered sphere must fit above the central ridge')


def volume_quadrature(scene, patch, lateral=10, vertical=8, sphere_radial=4,
                      sphere_azimuthal=16, sphere_vertical=2):
    """Integrate actual rectangular pits and the same sliced cylinders as S4."""
    patch.validate(scene)
    for value in (lateral, vertical, sphere_radial, sphere_azimuthal, sphere_vertical):
        if isinstance(value, bool) or not isinstance(value, int) or value < 2:
            raise ValueError('Source quadrature orders must be integers >= 2')
    def interval(n, lo, hi):
        x, w = np.polynomial.legendre.leggauss(n)
        return (lo+hi)/2+(hi-lo)*x/2, w*(hi-lo)/2
    points, weights, particle = [], [], []
    width = scene.pitch_um-scene.ridge_width_um
    x, wx = interval(lateral, -patch.length_um/2, patch.length_um/2)
    z, wz = interval(vertical, -scene.depth_um, 0)
    for center in patch.centers(scene):
        y, wy = interval(lateral, center-width/2, center+width/2)
        xx, yy, zz = np.meshgrid(x, y, z, indexing='ij')
        q = np.stack((xx, yy, zz), axis=-1).reshape(-1, 3)
        w = (wx[:, None, None]*wy[None, :, None]*wz[None, None, :]).ravel()
        points.append(q); weights.append(w); particle.append(np.zeros(len(w), bool))
    radius, dz = scene.diameter_um/2, scene.diameter_um/scene.slices
    phi = (np.arange(sphere_azimuthal)+.5)*2*np.pi/sphere_azimuthal
    for i in range(scene.slices):
        r = np.sqrt(radius**2-((i+.5)*dz-radius)**2)
        rho, wr = interval(sphere_radial, 0, r)
        # S4 adds sphere layers from top down, with midpoint disk radii.
        zz, wz = interval(sphere_vertical, scene.cap_top_um-(i+1)*dz, scene.cap_top_um-i*dz)
        rr, pp, zzz = np.meshgrid(rho, phi, zz, indexing='ij')
        q = np.stack((rr*np.cos(pp), rr*np.sin(pp), zzz), axis=-1).reshape(-1, 3)
        w = np.broadcast_to((rho*wr)[:, None, None]*wz[None, None, :]
                            *2*np.pi/sphere_azimuthal, rr.shape).ravel().copy()
        points.append(q); weights.append(w); particle.append(np.ones(len(w), bool))
    return np.concatenate(points), np.concatenate(weights), np.concatenate(particle)
