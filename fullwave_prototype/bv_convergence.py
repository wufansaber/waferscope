"""SCATMECH-side multipole convergence; run under NIST2, not the FEM interpreter."""
import argparse
import json
from pathlib import Path

import numpy as np


def main():
    from pySCATMECH.local import Local_BRDF_Model
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', default='build/theory_audit/bv_convergence')
    args = parser.parse_args()
    output = Path(args.output)
    output.mkdir(parents=True, exist_ok=True)
    theta = np.deg2rad(np.array((1, 15, 30, 45, 60, 75, 85)))
    phi = np.arange(16)*2*np.pi/16
    results = []
    for diameter, wavelength in ((100, 355), (30, 266), (200, 266)):
        amplitudes = []
        for order in (10, 14, 18):
            model = Local_BRDF_Model('Bobbert_Vlieger_BRDF_Model',
                {'lambda': wavelength/1000, 'radius': diameter/2000, 'sphere': complex(1.59, 0),
                 'substrate': complex(5.57, .387), 'delta': .01, 'type': 0,
                 'lmax': order, 'order': -1, 'Norm_Inc_Approx': 0})
            amplitudes.append(np.array([model.JonesDSC(np.deg2rad(75), t, f, coords='psps')
                                       for t in theta for f in phi]))
        fine = amplitudes[-1]
        results.append(dict(diameter_nm=diameter, wavelength_nm=wavelength,
            lmax10_vs18_complex_relative_l2=float(np.linalg.norm(amplitudes[0]-fine)/np.linalg.norm(fine)),
            lmax14_vs18_complex_relative_l2=float(np.linalg.norm(amplitudes[1]-fine)/np.linalg.norm(fine)),
            note='Fixed optical constants for numerical convergence, not measured dispersion. Both Jones input columns.'))
    report = dict(cases=results, angles=112, lmax=(10, 14, 18), gap_nm=10)
    (output/'report.json').write_text(json.dumps(report, indent=2), encoding='utf-8')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
