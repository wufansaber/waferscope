"""Combine completed cases from guarded runs without hiding an interrupted solver."""
import argparse
import json
from pathlib import Path

from fullwave_prototype.validation_plotting import plot_local
from fullwave_prototype.hybrid_run import load_sources


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--inputs', nargs='+', required=True)
    parser.add_argument('--guards', nargs='+', required=True)
    parser.add_argument('--output', required=True)
    args = parser.parse_args()
    reports = [json.loads((Path(p)/'validation_report.json').read_text(encoding='utf-8')) for p in args.inputs]
    mother = reports[0]['mother_sha256']
    if any(r['mother_sha256'] != mother or r['source'] != reports[0]['source'] for r in reports):
        raise ValueError('Do not combine different surfaces or physical sources')
    cases = [case for report in reports for case in report['cases']]
    _, metadata = load_sources(reports[0]['source'])
    rms = metadata['stationary_roughness']['rms_nm']
    for case in cases:
        case.setdefault('input_stationary_rms_nm', rms*case['height_scale'])
    result = dict(model='shared_surface_low_height_fem_spm_validation',
        status='complete_with_physical_limits', mother_sha256=mother, source=reports[0]['source'],
        cases=cases, full_model_converged=False,
        input_stationary_rms_nm=rms,
        resources=[json.loads(Path(p).read_text(encoding='utf-8')) for p in args.guards],
        scope='Completed selected-direction haze-only cases. One P2 UMFPACK attempt was stopped by '
              'memory protection; P2 sparsecholesky was completed separately. '
              'P1 and P2 denote FEM polynomial order, NOT perturbation-theory order. '
              'No calibration gain is applied; NOT full-model certification.')
    output = Path(args.output)
    output.mkdir(parents=True, exist_ok=False)
    (output/'validation_report.json').write_text(json.dumps(result, indent=2), encoding='utf-8')
    plot_local(output, result)


if __name__ == '__main__':
    main()
