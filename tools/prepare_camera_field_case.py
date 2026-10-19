"""Create a reproducible camera-matched 100 nm / 266 nm hybrid request."""
import argparse
from dataclasses import asdict, replace
import json
from pathlib import Path

from fullwave_prototype.imaging_run import build_parser, scattering_settings
from fullwave_prototype.imaging_options import imaging_settings
from fullwave_prototype.field_coverage import illumination_coverage


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', default='build/field_coverage_20261016/input.json')
    args = parser.parse_args()
    physical = build_parser().parse_args([
        '--wavelength-nm', '266', '--diameter-nm', '100', '--beam-profile', 'gaussian',
        '--beam-definition', 'wafer_gaussian', '--beam-order', '31',
        '--roughness-model', 'gaussian', '--rough-grid-points', '9', '--order', '1',
        '--accuracy', 'preview', '--mesh-scale', '1.5', '--pml-inner-fraction', '0',
        '--direct-solver', 'umfpack', '--threads', '2', '--illumination-fit', 'camera',
        '--minimum-field-illumination', '.5', '--illumination-power-w', '1',
        '--pupil-samples', '128', '--padding', '8', '--detector',
        '--detector-exposure-ms', '.01', '--detector-prnu-fraction', '.01', '--detector-dsnu-e', '.5',
        '--transmission', '.0001'])
    p, roughness = scattering_settings(physical)
    p = replace(p, domain_radius_um=.22, pml_thickness_um=.16, max_dofs=350000)
    settings = imaging_settings(physical, parameters=p)
    inputs = dict(parameters=asdict(p), roughness=asdict(roughness), imaging_settings=asdict(settings),
        coverage=.99, variance_fraction=.9, taper_width_um=.025, quadrature_order=5, mother_range_factor=2,
        assumptions=['Material n,k and camera QE remain assumed, not measured at 266 nm.',
                     'Optical transmission 1e-4 is explicit attenuation, not reduced source power.',
                     'One coherent realization; local FEM order 1 plus whole SPM, not certified full-wave accuracy.'])
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(inputs, indent=2), encoding='utf-8')
    output.with_name('imaging_settings.json').write_text(json.dumps(asdict(settings), indent=2), encoding='utf-8')
    check = illumination_coverage(p, settings.detector, settings.magnification,
        settings.numerical_apertures, settings.tilt_degrees, settings.minimum_field_illumination)
    print(json.dumps(dict(input=str(output.resolve()), length_y_um=p.beam_length_um,
        width_x_um=p.beam_width_um, irradiance_w_um2=settings.irradiance_w_um2, coverage=check), indent=2))


if __name__ == '__main__':
    main()
