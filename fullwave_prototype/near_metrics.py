"""Interface-aware diagnostics; sampled slice peaks are not 3-D field maxima."""
import numpy as np


def surface_distance(points, triangles):
    """Distance to a triangulated square height patch plus the exterior plane."""
    points, triangles = np.asarray(points, float), np.asarray(triangles, float)
    low, high = triangles[..., :2].min(axis=(0, 1)), triangles[..., :2].max(axis=(0, 1))
    edge = np.maximum(0, np.minimum(points[:, :2]-low, high-points[:, :2]).min(axis=1))
    best = points[:, 2]**2+edge**2
    for a, b, c in triangles:
        ab, ac = b-a, c-a
        normal = np.cross(ab, ac)
        normal /= np.linalg.norm(normal)
        delta = points-a
        height = delta@normal
        projected = delta-height[:, None]*normal
        aa, bb, cc = ab@ab, ab@ac, ac@ac
        u = (cc*(projected@ab)-bb*(projected@ac))/(aa*cc-bb*bb)
        v = (aa*(projected@ac)-bb*(projected@ab))/(aa*cc-bb*bb)
        inside = (u >= -1e-12) & (v >= -1e-12) & (u+v <= 1+1e-12)
        best = np.minimum(best, np.where(inside, height**2, np.inf))
        for start, end in ((a, b), (b, c), (c, a)):
            tangent = end-start
            fraction = np.clip((points-start)@tangent/(tangent@tangent), 0, 1)
            displacement = points-start-fraction[:, None]*tangent
            best = np.minimum(best, np.sum(displacement**2, axis=1))
    return np.sqrt(np.maximum(best, 0))


def clearances(points, p, triangles):
    radius = p.diameter_nm/2000
    center = np.array((0, 0, radius+p.gap_nm/1000))
    return surface_distance(points, triangles), abs(np.linalg.norm(points-center, axis=1)-radius)


def peak(points, electric, mask=None):
    intensity = np.sum(abs(electric)**2, axis=-1)
    mask = np.ones(len(points), bool) if mask is None else np.asarray(mask, bool)
    if not mask.any():
        raise ValueError('No samples remain in the requested region')
    indices = np.flatnonzero(mask)
    index = indices[np.argmax(intensity[mask])]
    return dict(intensity=float(intensity[index]), point_nm=(points[index]*1000).tolist(),
                retained_points=int(mask.sum()))


def comparison(points, base, refined, mask):
    if not np.any(mask):
        raise ValueError('Cannot compare an empty near-field region')
    a, b = base[mask], refined[mask]
    ia, ib = np.sum(abs(a)**2, axis=-1), np.sum(abs(b)**2, axis=-1)
    norm_e, norm_i = np.linalg.norm(b), np.linalg.norm(ib)
    if min(norm_e, norm_i) == 0:
        raise ValueError('Relative comparison requires a nonzero reference field')
    pa, pb = peak(points, base, mask), peak(points, refined, mask)
    return dict(complex_e_relative_l2=float(np.linalg.norm(a-b)/norm_e),
                intensity_relative_l2=float(np.linalg.norm(ia-ib)/norm_i),
                sampled_peak_relative_change=abs(pa['intensity']/pb['intensity']-1),
                base_peak=pa, refined_peak=pb)
