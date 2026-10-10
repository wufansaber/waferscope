"""Prepare small end-to-end haze checks for the new finite-envelope definitions."""
import argparse
from dataclasses import asdict, replace
import json
from pathlib import Path
from fullwave_prototype.parameters import Parameters
from fullwave_prototype.beam_inputs import fingerprint, wafer_dimensions
from fullwave_prototype.hybrid_run import optical_settings


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', required=True)
    args = parser.parse_args()
    root = Path(args.root)
    inputs = json.loads((root/'large_hybrid/hybrid_input.json').read_text(encoding='utf-8'))
    path = str((root/'beam_inputs/synthetic_intensity_demo.npz').resolve())
    for definition in ('wafer_rectangle', 'measured_wafer'):
        p = replace(Parameters(**inputs['parameters']), beam_definition=definition, beam_order=61,
            beam_measurement_path=path if definition == 'measured_wafer' else '',
            beam_measurement_sha256=fingerprint(path) if definition == 'measured_wafer' else '')
        length, width = wafer_dimensions(p)
        p = replace(p, beam_length_um=length, beam_width_um=width)
        values = dict(inputs['imaging_settings'], illumination_definition=definition,
                      illumination_length_um=length, illumination_width_um=width)
        settings = optical_settings(values, p)
        output = root/('haze_'+definition+'_input.json')
        output.write_text(json.dumps(dict(inputs, parameters=asdict(p), imaging_settings=asdict(settings),
            samples=1, target_relative_ci=.05), indent=2), encoding='utf-8')
        print(output)


if __name__ == '__main__':
    main()
