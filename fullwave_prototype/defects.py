"""Finite raised defects, shared by CAD, input validation and field figures."""
import math
from types import SimpleNamespace

import numpy as np


DEFAULTS = dict(defect_shape='sphere', defect_length_nm=200., defect_width_nm=50.,
                defect_height_nm=50., defect_rotation_deg=0.)
SHAPES = ('sphere', 'line_x', 'line_y', 'cross')


def values(p):
    return SimpleNamespace(**dict(DEFAULTS, **p)) if isinstance(p, dict) else p


def footprint_nm(p):
    p = values(p)
    if p.defect_shape == 'sphere':
        raise ValueError('A sphere has no prismatic footprint')
    a, b = p.defect_length_nm/2, p.defect_width_nm/2
    if p.defect_shape == 'cross' and a > b:
        points = [(-b, -a), (b, -a), (b, -b), (a, -b), (a, b), (b, b),
                  (b, a), (-b, a), (-b, b), (-a, b), (-a, -b), (-b, -b)]
    else:
        points = [(-a, -b), (a, -b), (a, b), (-a, b)]
    angle = math.radians(p.defect_rotation_deg+(90 if p.defect_shape == 'line_y' else 0))
    c, s = math.cos(angle), math.sin(angle)
    return np.asarray(points) @ np.array([[c, s], [-s, c]])


def center_height_um(p):
    p = values(p)
    height = p.diameter_nm if p.defect_shape == 'sphere' else p.defect_height_nm
    return (p.gap_nm+height/2)/1000


def required_radius_um(p):
    p = values(p)
    if p.defect_shape == 'sphere':
        return (p.diameter_nm+p.gap_nm)/1000
    return float(np.sqrt(np.max(np.sum(footprint_nm(p)**2, axis=1))+
                         (p.gap_nm+p.defect_height_nm)**2)/1000)


def validate_geometry(p):
    if p.defect_shape not in SHAPES:
        raise ValueError('Defect shape must be sphere, line_x, line_y or cross')
    if min(p.defect_length_nm, p.defect_width_nm, p.defect_height_nm) <= 0:
        raise ValueError('Defect length, width and height must be positive')
    if p.defect_shape != 'sphere' and p.defect_width_nm > p.defect_length_nm:
        raise ValueError('Defect width cannot exceed its length')
    minimum = required_radius_um(p)
    if p.domain_radius_um <= minimum:
        raise ValueError('The defect must be entirely inside the physical domain; '
                         'domain radius must exceed {:.6g} um'.format(minimum))


def description(p):
    p = values(p)
    if p.defect_shape == 'sphere':
        return 'sphere {:g} nm'.format(p.diameter_nm)
    return '{} {:g} x {:g} x {:g} nm, rotation {:g} deg'.format(
        p.defect_shape, p.defect_length_nm, p.defect_width_nm,
        p.defect_height_nm, p.defect_rotation_deg)


def metadata(p):
    p = values(p)
    if p.defect_shape == 'sphere':
        volume = 4*np.pi/3*(p.diameter_nm/2000)**3
    else:
        area = p.defect_length_nm*p.defect_width_nm
        if p.defect_shape == 'cross':
            area = 2*area-p.defect_width_nm**2
        volume = area*p.defect_height_nm/1e9
    return dict(shape=p.defect_shape, description=description(p),
                expected_volume_um3=float(volume), bottom_z_nm=p.gap_nm,
                center_height_um=center_height_um(p), required_radius_um=required_radius_um(p),
                scope='Finite solid above nominal wafer; not an etched groove or infinite grating')


def cad_solid(p, scale=1.):
    from netgen.occ import Sphere, Face, MakePolygon, Vertex, Vec
    if p.defect_shape == 'sphere':
        return Sphere((0, 0, center_height_um(p)*scale), p.diameter_nm*scale/2000)
    xy = footprint_nm(p)*scale/1000
    vertices = [Vertex((x, y, p.gap_nm*scale/1000)) for x, y in xy]
    # One concave outline: the cross intersection is material only once.
    face = Face(MakePolygon(vertices+[vertices[0]]))
    return face.Extrude(Vec(0, 0, p.defect_height_nm*scale/1000))


def boundary_distance_um(points, p):
    p = values(p)
    points = np.asarray(points, float)
    if p.defect_shape == 'sphere':
        return abs(np.linalg.norm(points-[0, 0, center_height_um(p)], axis=-1)-p.diameter_nm/2000)
    polygon = footprint_nm(p)/1000
    xy = points[..., :2]
    distance2 = np.full(points.shape[:-1], np.inf)
    inside = np.zeros(points.shape[:-1], bool)
    for a, b in zip(polygon, np.roll(polygon, -1, axis=0)):
        edge = b-a
        t = np.clip(np.sum((xy-a)*edge, axis=-1)/(edge@edge), 0, 1)
        distance2 = np.minimum(distance2, np.sum((xy-a-t[..., None]*edge)**2, axis=-1))
        if edge[1] != 0:
            inside ^= ((a[1] > xy[..., 1]) != (b[1] > xy[..., 1])) & (
                xy[..., 0] < a[0]+(xy[..., 1]-a[1])*edge[0]/edge[1])
    lateral = np.sqrt(distance2)*np.where(inside, -1, 1)
    vertical = abs(points[..., 2]-center_height_um(p))-p.defect_height_nm/2000
    q = np.stack((lateral, vertical), axis=-1)
    signed = np.linalg.norm(np.maximum(q, 0), axis=-1)+np.minimum(np.max(q, axis=-1), 0)
    return abs(signed)


def section_nm(p, plane, offset_nm=0.):
    """Prism outline in XY or a central vertical slice, in plotted coordinates."""
    p = values(p)
    polygon = footprint_nm(p)
    if plane == 'xy':
        return polygon if p.gap_nm <= offset_nm <= p.gap_nm+p.defect_height_nm else None
    fixed, horizontal = (1, 0) if plane == 'xz' else (0, 1)
    hits = []
    for a, b in zip(polygon, np.roll(polygon, -1, axis=0)):
        da, db = a[fixed]-offset_nm, b[fixed]-offset_nm
        if abs(da) < 1e-10:
            hits.append(a[horizontal])
        if da*db < 0:
            hits.append((a+(b-a)*da/(da-db))[horizontal])
    if len(hits) < 2:
        return None
    low, high = min(hits), max(hits)
    bottom, top = p.gap_nm, p.gap_nm+p.defect_height_nm
    return np.array([(low, bottom), (high, bottom), (high, top), (low, top)])


def draw_section(ax, p, plane, offset_nm=0., **style):
    from matplotlib.patches import Circle, Polygon
    p = values(p)
    style = dict(dict(fill=False, edgecolor='cyan', linewidth=1), **style)
    if p.defect_shape == 'sphere':
        center = np.array([0., 0., center_height_um(p)*1000])
        horizontal, vertical, fixed = {'xz': (0, 2, 1), 'xy': (0, 1, 2), 'yz': (1, 2, 0)}[plane]
        radius2 = (p.diameter_nm/2)**2-(offset_nm-center[fixed])**2
        if radius2 > 0:
            ax.add_patch(Circle(center[[horizontal, vertical]], np.sqrt(radius2), **style))
    else:
        outline = section_nm(p, plane, offset_nm)
        if outline is not None:
            ax.add_patch(Polygon(outline, closed=True, **style))
