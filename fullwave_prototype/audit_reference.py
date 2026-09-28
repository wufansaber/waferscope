"""Independent NIST2 reference for the native Maxwell validation harness."""
import argparse
import json
from pathlib import Path
import numpy as np


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--request', required=True)
    parser.add_argument('--output', required=True)
    args = parser.parse_args()
    request = json.loads(Path(args.request).read_text(encoding='utf-8'))
    p = request['parameters']
    with np.load(request['directions_file']) as data:
        directions, weights = data['directions'], data['weights']
    if request['kind'] == 'mie':
        from pySCATMECH.scatterer import Free_Space_Scatterer
        model = Free_Space_Scatterer('MieScatterer', {'lambda': p['wavelength_nm']/1000,
            'radius': p['diameter_nm']/2000, 'sphere': complex(p['particle_n'], p['particle_k']), 'medium': 1})
        t = np.deg2rad(p['incidence_deg'])
        incident = (0, np.sin(t), -np.cos(t))
        dsc = np.array([np.asarray(model.DifferentialScatteringCrossSection(incident, k))[0, 0]
                        for k in directions])
        np.savez_compressed(args.output, csca_um2=np.sum(weights*dsc),
                            cext_um2=np.asarray(model.Extinction(incident))[0, 0])
    else:
        from pySCATMECH.local import Local_BRDF_Model
        model = Local_BRDF_Model('Bobbert_Vlieger_BRDF_Model', {'lambda': p['wavelength_nm']/1000,
            'radius': p['diameter_nm']/2000, 'sphere': complex(p['particle_n'], p['particle_k']),
            'substrate': complex(p['substrate_n'], p['substrate_k']), 'delta': p['gap_nm']/1000,
            'type': 0, 'lmax': 14, 'order': -1, 'Norm_Inc_Approx': 0})
        column = 0 if p['polarization'] == 's' else 1
        amplitudes = np.array([model.JonesDSC(np.deg2rad(p['incidence_deg']),
            np.arccos(np.clip(k[2], -1, 1)), np.arctan2(k[1], k[0])-np.pi/2,
            coords='psps')[:, column] for k in directions])
        np.savez_compressed(args.output, amplitude=amplitudes)


if __name__ == '__main__':
    main()
