"""Manual-sourced Gemi camera assumptions for a single-exposure area approximation."""
import argparse
from dataclasses import asdict
import json
from pathlib import Path

from detector_model import DetectorParameters
from camera_profiles import gemi_area_parameters


MANUAL = 'D:/项目/相机/Gemi 8KTDI 用户手册V1.0.0_2025032600.pdf'


def profile(gain_e_dn=16., dark_current_e_s=1000.):
    p = gemi_area_parameters(gain_e_dn, dark_current_e_s)
    p.validate()
    return dict(detector=asdict(p), provenance=dict(manual=MANUAL,
        appendix='Appendix 3, printed pages 67-68 / PDF pages 68-69',
        mode='Single area exposure approximation, NO TDI stage multiplication',
        manual_values=dict(pixel_pitch_um=5., qe_266nm_lower_bound=.639,
            read_noise_median_e_10bit=12.1, full_well_lower_bound_e=16000,
            dark_current_upper_e_s_at_10C=1000, dark_current_upper_e_s_at_40C=7000,
            corrected_dsnu_upper_e=5.6, uncorrected_dsnu_upper_e=10.8, prnu_upper_fraction=.00124),
        assumptions=['The manual states Area mode is for focusing and image quality is not guaranteed (printed page 27).',
            '10-bit noise specifications are adopted as hypothetical area-mode parameters, not guaranteed area measurements.',
            'Conversion gain e-/DN is NOT given. 16 e-/DN is an estimate near 16000/1023, not a calibration.',
            'Black level 64 DN is an assumed operating setting. Analog/digital gain mapping is not calibrated.',
            '12.1 e- median is used as uniform per-pixel Gaussian input read-noise sigma; ADC rounding is also simulated. '
            'The manual does not specify whether its reported noise already contains ADC quantization.',
            'QE, full well, dark current, DSNU and PRNU bounds are used as numerical inputs; they are not exact measurements.',
            '20x20 pixels is a 100x100 um computational sensor ROI, not a claim of legal hardware ROI settings.',
            'Exposure 10 us; paired same-surface reference; fixed patterns persist across frames.']))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--gain-e-dn', type=float, default=16.)
    parser.add_argument('--dark-current-e-s', type=float, default=1000.)
    parser.add_argument('--output', default='build/field_coverage_20261016/gemi_area/camera_profile.json')
    args = parser.parse_args()
    values = profile(args.gain_e_dn, args.dark_current_e_s)
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(values, indent=2), encoding='utf-8')
    print(str(output.resolve()))


if __name__ == '__main__':
    main()
