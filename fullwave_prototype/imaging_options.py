"""Shared imaging-only CLI options, also used when reusing scattering sources."""
from detector_model import DetectorParameters
from .imaging import ImagingSettings


IMAGING_FLAGS = {'na', 'tilts', 'magnification', 'objective-focal-mm', 'irradiance-w-um2',
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


def imaging_settings(args):
    detector = DetectorParameters(enabled=args.detector, **{
        key.replace('-', '_'): getattr(args, 'detector_'+key.replace('-', '_'))
        for key in DETECTOR_OPTIONS})
    settings = ImagingSettings(numerical_apertures=tuple(args.na), tilt_degrees=tuple(args.tilts),
        magnification=args.magnification, focal_length_mm=args.objective_focal_mm,
        irradiance_w_um2=args.irradiance_w_um2, image_defocus_um=args.defocus_um,
        pupil_samples=args.pupil_samples, padding=args.padding,
        transmission=args.transmission, analyzer=args.analyzer,
        aberrations_nm=(args.astigmatism_nm, args.coma_nm, args.spherical_nm), detector=detector)
    settings.validate()
    return settings
