"""Run with NIST2 to produce an independent SCATMECH reference, without installing FEM there."""
import argparse
import json
from pathlib import Path
import numpy as np

from .parameters import Parameters
from .far_field import angular_grid


def calculate(p, theta, phi):
    if p.defect_shape != 'sphere':
        raise ValueError('Bobbert-Vlieger reference supports spheres only, not line or cross defects')
    from pySCATMECH.local import Local_BRDF_Model
    model = Local_BRDF_Model('Bobbert_Vlieger_BRDF_Model', {
        'lambda': p.wavelength_nm/1000, 'radius': p.diameter_nm/2000,
        'sphere': complex(p.particle_n, p.particle_k),
        'substrate': complex(p.substrate_n, p.substrate_k),
        'delta': p.gap_nm/1000, 'type': 0, 'lmax': 14,
        'order': -1, 'Norm_Inc_Approx': 0})
    index = 0 if p.polarization == 's' else 1
    amplitudes = np.zeros((len(theta), len(phi), 2), dtype=complex)
    for i, t in enumerate(theta):
        for j, f in enumerate(phi):
            amplitudes[i, j] = np.asarray(model.JonesDSC(
                np.deg2rad(p.incidence_deg), min(t, np.pi/2-1e-8), f-np.pi/2,
                coords='psps'))[:, index]
    return amplitudes


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--config', required=True)
    parser.add_argument('--output', required=True)
    args = parser.parse_args()
    p = Parameters(**json.loads(Path(args.config).read_text())['parameters'])
    theta, phi, _ = angular_grid(5)
    amplitudes = calculate(p, theta, phi)
    np.savez_compressed(args.output, theta=theta, phi=phi, amplitude_sp_um=amplitudes,
                        dsc_um2_sr=np.sum(abs(amplitudes)**2, axis=-1))
