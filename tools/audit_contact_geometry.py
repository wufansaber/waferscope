"""Record the actual SCATMECH separation and FEM contact geometry for a release."""
import argparse
from dataclasses import replace
import hashlib
import json
from pathlib import Path

from fullwave_prototype.defects import metadata
from fullwave_prototype.imaging_run import build_parser, scattering_settings
from scatmech_ui.parameters import SimulationParameters
from scatmech_ui.surface_model import build_scattering_model


def main():
    from pySCATMECH.local import Local_BRDF_Model
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    p, roughness = scattering_settings(build_parser().parse_args([]))
    default_delta = float(Local_BRDF_Model('Bobbert_Vlieger_BRDF_Model').getParameters()['delta'])
    rows = []
    for diameter in range(10, 101, 10):
        params = build_scattering_model(SimulationParameters(particle_diameter_nm=diameter)).model.getParameters()
        gap_um, radius_um = float(params['delta']), float(params['radius'])
        geometry = metadata(replace(p, diameter_nm=diameter))
        if gap_um != 0 or geometry['bottom_z_nm'] != 0 or abs(radius_um-diameter/2000) > 1e-14:
            raise AssertionError('The sphere is not in contact')
        rows.append(dict(diameter_nm=diameter, scatmech_delta_um=gap_um, scatmech_radius_um=radius_um,
                         fem_geometry=geometry))
    root = Path(__file__).resolve().parents[1]
    files = ('scatmech_ui/surface_model.py', 'fullwave_prototype/imaging_run.py',
             'fullwave_prototype/rough_run.py', 'fullwave_ui/controls.py', 'fullwave_ui/defect_controls.py',
             'scatmech_ui/imaging.py', 'scatmech_ui/field_reference.py', 'fullwave_prototype/solver.py')
    report = dict(status='passed', library_default_delta_um=default_delta, rows=rows,
        fem_default_gap_nm=p.gap_nm, fem_default_rms_nm=roughness.rms_nm,
        propagated_field_reference='top-wafer plane z=0', native_jones_reference='sphere center',
        preserved_rough_preset=dict(gap_nm=10, rms_nm=.2, contact=False),
        source_sha256={name: hashlib.sha256((root/name).read_bytes()).hexdigest() for name in files},
        scope='Geometry and defaults only. Not a physics-accuracy certification. Old archives are not modified.')
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2), encoding='utf-8')
    print('Contact geometry: {} sizes passed; library delta={:g} um; FEM gap={:g} nm'.format(
        len(rows), default_delta, p.gap_nm))


if __name__ == '__main__':
    main()
