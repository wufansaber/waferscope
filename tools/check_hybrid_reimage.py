"""Check real exported hybrid power scaling and disabled apertures."""
import argparse
import json
from pathlib import Path

import numpy as np


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--baseline', required=True)
    parser.add_argument('--reimage', required=True)
    args = parser.parse_args()
    baseline, reimage = Path(args.baseline), Path(args.reimage)
    reports = [json.loads((directory/'hybrid_report.json').read_text(encoding='utf-8'))
               for directory in (baseline, reimage)]
    assert all(r['model'] == 'whole_spm_local_fem_hybrid' for r in reports)
    assert reports[0]['mother_surface_sha256'] == reports[1]['mother_surface_sha256']
    scale = reports[1]['settings']['illumination_power_w']/reports[0]['settings']['illumination_power_w']
    with np.load(baseline/'image_fields.npz', allow_pickle=False) as original, np.load(
            reimage/'image_fields.npz', allow_pickle=False) as changed:
        for index, (before, after) in enumerate(zip(reports[0]['channels'], reports[1]['channels']), 1):
            for state in ('with', 'without', 'change', 'interference'):
                key = 'ch{}_image_{}_w_um2'.format(index, state)
                if after['na'] == 0:
                    np.testing.assert_array_equal(changed[key], 0)
                else:
                    assert before['na'] == after['na'] and before['tilt_deg'] == after['tilt_deg']
                    np.testing.assert_allclose(changed[key], scale*original[key], rtol=1e-12, atol=1e-25)
            if after['na'] == 0:
                assert after['with']['image_power_w'] == 0
            else:
                assert np.isclose(after['with']['image_power_w'], scale*before['with']['image_power_w'], rtol=1e-12)
            assert after['decomposition_relative_max_error'] < 1e-12
    print('Same mother preserved; absolute power scaling, NA=0 and coherent decomposition passed')


if __name__ == '__main__':
    main()
