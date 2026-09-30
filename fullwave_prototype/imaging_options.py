"""Shared imaging-only CLI options, also used when reusing scattering sources."""
import numpy as np

from detector_model import DetectorParameters
from .imaging import ImagingSettings
from .parameters import Parameters


IMAGING_FLAGS = {'na', 'tilts', 'magnification', 'objective-focal-mm', 'irradiance-w-um2',
                 'illumination-power-w', 'illumination-length-um', 'illumination-width-um',
                 'defocus-um', 'transmission', 'analyzer', 'pupil-samples', 'padding',
                 'archive-mode', 'astigmatism-nm', 'coma-nm', 'spherical-nm', 'detector'}
DETECTOR_OPTIONS = {'pixel-pitch-um': 6.5, 'exposure-ms': 1., 'quantum-efficiency': .6,
                    'read-noise-e': 2., 'dark-current-e-s': .1, 'background-e-s': 0.,
                    'full-well-e': 30000., 'pixels-x': 16, 'pixels-y': 16, 'seed': 0}
IMAGING_FLAGS.update('detector-'+key for key in DETECTOR_OPTIONS)


def add_options(parser):
    parser.add_argument('--archive-mode', choices=('compact', 'full'), default='compact')
    parser.add_argument('--pupil-cache', help='Optional content-addressed complex pupil cache directory')
    for key in ('astigmatism', 'coma', 'spherical'):
        parser.add_argument('--'+key+'-nm', type=float, default=0)
    parser.add_argument('--detector', action='store_true')
    for key, value in DETECTOR_OPTIONS.items():
        parser.add_argument('--detector-'+key, type=type(value), default=value)


def imaging_settings(args, incidence_deg=None):
    detector = DetectorParameters(enabled=args.detector, **{
        key.replace('-', '_'): getattr(args, 'detector_'+key.replace('-', '_'))
        for key in DETECTOR_OPTIONS})
    footprint = (args.illumination_power_w, args.illumination_length_um,
                 args.illumination_width_um)
    if all(value is None for value in footprint) and args.irradiance_w_um2 is None:
        footprint = (1., 2., 2.)
    if any(value is not None for value in footprint):
        if any(value is None for value in footprint):
            raise ValueError('Specify illumination power, length and width together')
        incidence = (incidence_deg if incidence_deg is not None else
                     args.incidence_deg if args.incidence_deg is not None else
                     Parameters().incidence_deg)
        if incidence is None or not np.isfinite(incidence) or not 0 <= incidence < 90:
            raise ValueError('Source incidence is required for a wafer-surface illumination footprint')
        power, length, width = footprint
        if not all(np.isfinite(value) and value > 0 for value in footprint):
            raise ValueError('Illumination power and footprint dimensions must be positive and finite')
        irradiance = power/(length*width*np.cos(np.deg2rad(incidence)))
        if args.irradiance_w_um2 is not None and not np.isclose(
                args.irradiance_w_um2, irradiance, rtol=1e-10, atol=0):
            raise ValueError('Explicit irradiance conflicts with illumination power and footprint')
    else:
        power = length = width = incidence = None
        irradiance = args.irradiance_w_um2 if args.irradiance_w_um2 is not None else 1e-6
    settings = ImagingSettings(numerical_apertures=tuple(args.na), tilt_degrees=tuple(args.tilts),
        magnification=args.magnification, focal_length_mm=args.objective_focal_mm,
        irradiance_w_um2=irradiance, illumination_power_w=power,
        illumination_length_um=length, illumination_width_um=width,
        illumination_incidence_deg=incidence, image_defocus_um=args.defocus_um,
        pupil_samples=args.pupil_samples, padding=args.padding,
        transmission=args.transmission, analyzer=args.analyzer,
        aberrations_nm=(args.astigmatism_nm, args.coma_nm, args.spherical_nm), detector=detector)
    settings.validate()
    return settings
