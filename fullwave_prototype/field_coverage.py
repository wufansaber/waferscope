"""Wafer footprints of independent cameras and unchanged-beam coverage checks."""
from dataclasses import replace

import numpy as np


FIT_MODES = ('manual', 'camera')


def camera_footprints(detector, magnification, apertures, tilts):
    detector.validate()
    if not np.isfinite(magnification) or magnification < 1:
        raise ValueError('Magnification must be finite and >= 1')
    if len(apertures) != 3 or len(tilts) != 3:
        raise ValueError('Three independent objective settings required')
    width = detector.pixels_x*detector.pixel_pitch_um/magnification
    height = detector.pixels_y*detector.pixel_pitch_um/magnification
    footprints = []
    for i, (na, tilt) in enumerate(zip(apertures, tilts)):
        if not np.isfinite((na, tilt)).all() or not 0 <= na < 1 or abs(tilt) > 89.9:
            raise ValueError('Invalid objective aperture or tilt')
        if na > 0:
            footprints.append(dict(channel=i+1, tilt_deg=float(tilt),
                width_x_um=float(width/np.cos(np.deg2rad(tilt))), height_y_um=float(height),
                wafer_edge_object_defocus_um=float(width/2*abs(np.tan(np.deg2rad(tilt))))))
    return footprints


def nominal_dimensions(footprints, minimum_relative_flux=.5):
    if not footprints:
        raise ValueError('Enable at least one objective to match illumination to a camera')
    if not np.isfinite(minimum_relative_flux) or not .05 <= minimum_relative_flux <= .95:
        raise ValueError('Minimum illumination / center must be in [0.05, 0.95]')
    # I/I0 = exp[-8(x/Dx)^2-8(y/Dy)^2]; both corner coordinates contribute.
    scale = np.sqrt(-4/np.log(minimum_relative_flux))
    return (float(max(f['height_y_um'] for f in footprints)*scale),
            float(max(f['width_x_um'] for f in footprints)*scale))


def _flux_map(p, footprint, order, samples=17):
    from .surface_beam import surface_fields_order
    x = np.linspace(-footprint['width_x_um']/2, footprint['width_x_um']/2, samples)
    y = np.linspace(-footprint['height_y_um']/2, footprint['height_y_um']/2, samples)
    _, e, h = surface_fields_order(x, y, p, order)
    return -np.real(np.cross(e, h.conj()))[..., 2]


def illumination_coverage(p, detector, magnification, apertures, tilts,
                          minimum_relative_flux=.5):
    footprints = camera_footprints(detector, magnification, apertures, tilts)
    channels = []
    for footprint in footprints:
        if p.beam_profile == 'gaussian':
            order = max(41, p.beam_order)
            coarse = _flux_map(p, footprint, order)
            fine = _flux_map(p, footprint, max(61, order+20))
            error = float(np.linalg.norm(coarse-fine)/max(np.linalg.norm(fine), 1e-300))
            center = float(fine[len(fine)//2, len(fine)//2])
            if not np.isfinite(center) or center <= 0:
                raise ValueError('Invalid incident wafer flux at the field center')
            ratio = float(fine.min()/center)
        elif p.beam_waist_um:
            ratio, error = None, None
        else:
            ratio, error = 1., 0.
        channels.append(dict(footprint, minimum_relative_incident_flux=ratio,
            field_quadrature_relative_l2=error,
            target_met=bool(ratio is not None and ratio >= minimum_relative_flux and error < .005)))
    return dict(channels=channels, target_minimum_relative_incident_flux=minimum_relative_flux,
        all_channels_covered=bool(channels and all(c['target_met'] for c in channels)),
        scope='Geometric sensor footprint projected to z=0; incident downward Poynting flux '
              'sampled on 17x17 points per footprint, not a rigorous continuous lower bound. '
              'Legacy beam-waist envelopes are not certified by this coverage check. '
              'Tilted wafer edges are defocused; coverage does not certify sharp imaging or uniform haze.')


def fit_camera_beam(p, detector, magnification, apertures, tilts, minimum_relative_flux=.5):
    if p.beam_profile != 'gaussian' or p.beam_definition not in ('wafer_gaussian', 'beam_cross_section'):
        raise ValueError('Camera matching requires a Gaussian beam; use manual for measured or rectangular beams')
    footprints = camera_footprints(detector, magnification, apertures, tilts)
    length, width = nominal_dimensions(footprints, minimum_relative_flux)
    if p.beam_definition == 'beam_cross_section':
        length *= np.cos(np.deg2rad(p.incidence_deg))
    candidate = replace(p, beam_length_um=length, beam_width_um=width)
    # Check the actual propagating vector field, not just its nominal Gaussian envelope.
    low, high = 1., 1.02
    report = illumination_coverage(candidate, detector, magnification, apertures, tilts, minimum_relative_flux)
    if not report['all_channels_covered']:
        for _ in range(12):
            candidate = replace(p, beam_length_um=length*high, beam_width_um=width*high)
            report = illumination_coverage(candidate, detector, magnification, apertures, tilts, minimum_relative_flux)
            if report['all_channels_covered']:
                break
            low, high = high, high*1.1
        else:
            raise ValueError('Camera illumination fit did not converge; inspect beam sampling and field size')
        while high-low > .002*high:
            mid = (low+high)/2
            trial = replace(p, beam_length_um=length*mid, beam_width_um=width*mid)
            check = illumination_coverage(trial, detector, magnification, apertures, tilts, minimum_relative_flux)
            if check['all_channels_covered']:
                high, candidate, report = mid, trial, check
            else:
                low = mid
    return candidate, report


def rough_surface_coverage(settings, plan):
    footprints = camera_footprints(settings.detector, settings.magnification,
                                  settings.numerical_apertures, settings.tilt_degrees)
    required = max((max(f['width_x_um'], f['height_y_um']) for f in footprints), default=0.)
    usable = plan['width_um']-2*plan['taper_width_um']
    report = dict(required_wafer_span_um=required, untapered_surface_span_um=usable,
                  camera_fields_inside_untapered_surface=bool(usable >= required))
    if settings.illumination_fit == 'camera' and usable < required:
        raise ValueError('Rough surface/taper does not cover the camera footprints; enlarge the surface')
    return report
