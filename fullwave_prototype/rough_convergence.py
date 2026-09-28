"""Compare saved rough-surface discretizations, not an exact roughness oracle."""
import json
from pathlib import Path
import numpy as np

from .parameters import Parameters
from .physical_audit import sphere_rule
from .far_field import radiation_amplitude
from .audit_optics import near_mesh_check


def main():
    root = Path('build/theory_audit')
    source = root/'rough_balanced/scattering'
    old = Path('build/boundary_upgrade/order3_condensed/rough_pml')
    def read(path):
        return json.loads(path.read_text(encoding='utf-8'))
    p = Parameters(**read(source/'report.json')['parameters'])
    old_parameters = read(old/'report.json')['parameters']
    if read(source/'report.json')['roughness'] != read(old/'report.json')['roughness']:
        raise ValueError('Compared roughness realizations differ')
    for key in ('wavelength_nm', 'diameter_nm', 'particle_n', 'particle_k', 'substrate_n',
                'substrate_k', 'incidence_deg', 'polarization', 'gap_nm', 'order'):
        if getattr(p, key) != old_parameters[key]:
            raise ValueError('Compared cases differ in '+key)
    directions, measure = sphere_rule(20, hemisphere=True)
    with np.load(source/'paired_fields.npz') as data:
        actual = radiation_amplitude(data['source_points_with'],
            data['source_polarization_weights_with'], directions, p)
    with np.load(old/'fields.npz') as data:
        reference = data['amplitude_1']
    weights = np.sqrt(measure)[:, None]
    actual_power = float(np.sum(measure*np.sum(abs(actual)**2, axis=1)))
    reference_power = float(np.sum(measure*np.sum(abs(reference)**2, axis=1)))
    result = dict(upper_dsc_um2=actual_power, reference_upper_dsc_um2=reference_power,
        relative_hemisphere_difference=abs(actual_power/reference_power-1),
        weighted_complex_amplitude_difference=float(np.linalg.norm((actual-reference)*weights)/np.linalg.norm(reference*weights)),
        near_order2_to_balanced=near_mesh_check('build/fullwave_imaging/complete_preview/scattering', source),
        reference=str(old.resolve()), reference_strength=1.,
        note='Order-3 original PML mesh vs balanced PML mesh at strength 1. Discretization sensitivity only, no independent rough-surface truth.')
    (root/'rough_convergence.json').write_text(json.dumps(result, indent=2, allow_nan=False), encoding='utf-8')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
