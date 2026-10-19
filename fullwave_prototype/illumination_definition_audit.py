"""Compare actual incident fields for explicitly different coherent beam inputs."""
import argparse
from dataclasses import replace, asdict
import json
from pathlib import Path
import numpy as np
from scipy.integrate import trapezoid
from .parameters import Parameters
from .beam_inputs import fingerprint, measured_envelope, wafer_dimensions, truncation_report
from .finite_illumination import power_area
from .surface_beam import surface_fields_order


def run(source, output, measurement=None, phase_assumption='flat', native_check=False):
    source, output = Path(source), Path(output)
    output.mkdir(parents=True, exist_ok=True)
    metadata = json.loads((source/'hybrid_source_report.json').read_text(encoding='utf-8'))
    p = replace(Parameters(**metadata['parameters']), beam_order=101)
    demo = measurement is None
    if demo:
        measurement = output/'synthetic_intensity_demo.npz'
        axis = (np.arange(32)-16)*.1875
        x, y = np.meshgrid(axis, axis)
        np.savez_compressed(measurement, x_um=axis, y_um=axis, intensity=np.exp(-2*(x*x+y*y)))
    path = str(Path(measurement).resolve())
    digest = fingerprint(path)
    _, _, _, info = measured_envelope(path, digest, phase_assumption)
    reports, images, axes = [], [], []
    if native_check:
        import ngsolve as ng
        from netgen.occ import Box, OCCGeometry
        mesh = ng.Mesh(OCCGeometry(Box((-.3, -.3, -.3), (.3, .3, .3))).GenerateMesh(maxh=.2))
    for definition in ('wafer_gaussian', 'beam_cross_section', 'wafer_rectangle', 'measured_wafer'):
        beam = replace(p, beam_definition=definition, beam_measurement_path=path if definition == 'measured_wafer' else '',
            beam_measurement_sha256=digest if definition == 'measured_wafer' else '', beam_phase_assumption=phase_assumption)
        if definition == 'measured_wafer':
            beam = replace(beam, beam_length_um=info['input_extent_y_um'], beam_width_um=info['input_extent_x_um'])
        beam.validate()
        length, width = wafer_dimensions(beam)
        x, y = np.linspace(-1.5*width, 1.5*width, 128), np.linspace(-1.5*length, 1.5*length, 128)
        _, e, h = surface_fields_order(x, y, beam, 101)
        checks, previous, stable, final_order = [], e, 0, 101
        for order in (151, 251, 401, 601):
            _, fine_e, fine_h = surface_fields_order(x, y, beam, order)
            field_change = float(np.linalg.norm(previous-fine_e)/max(np.linalg.norm(fine_e), 1e-300))
            checks.append(dict(orders=[final_order, order], relative_l2_change=field_change))
            stable = stable+1 if field_change <= .001 else 0
            previous, final_order = fine_e, order
            if stable >= 2:
                break
        irradiance = 1/power_area(beam)
        flux = -np.real(np.cross(fine_e, fine_h.conj()))[..., 2]*irradiance
        report = truncation_report(beam)
        report.update(parameters=asdict(beam), total_laser_power_w=1.,
            field_quadrature_orders=checks[-1]['orders'], field_quadrature_relative_change=field_change,
            field_quadrature_checks=checks, field_quadrature_target_met=stable >= 2,
            displayed_incident_power_w=float(trapezoid(trapezoid(flux, x, axis=1), y)),
            measurement_is_synthetic_demo=demo if definition == 'measured_wafer' else False)
        if native_check:
            from .solver import background_coefficient
            from .background_sampling import electric_background
            probe_beam = replace(beam, beam_order=21)
            # Normal E is discontinuous at z=0; compare strictly inside each medium.
            xyz = np.array([[.1, .2, .1], [.05, .02, -.1], [.1, -.1, .04]])
            locations = mesh(xyz[:, 0], xyz[:, 1], xyz[:, 2])
            reference = np.asarray(background_coefficient(probe_beam)(locations)).reshape(-1, 3)
            batch = electric_background(xyz, probe_beam)
            report['native_background_comparison_relative_l2'] = float(
                np.linalg.norm(batch-reference)/max(np.linalg.norm(reference), 1e-300))
        reports.append(report)
        images.append(flux)
        axes.append((x, y))
    from .imaging_plotting import _pyplot
    plt = _pyplot()
    fig, plots = plt.subplots(2, 2, figsize=(10, 9), layout='constrained')
    for ax, flux, (x, y), report in zip(plots.flat, images, axes, reports):
        image = ax.pcolormesh(x, y, flux, cmap='viridis', shading='auto')
        ax.set(xlabel='Wafer X (um)', ylabel='Wafer Y (um)', aspect='equal',
            title='{}\nExcluded scalar spectrum: {:.3%}'.format(report['definition'], report['excluded_scalar_spectral_fraction']))
        fig.colorbar(image, ax=ax, label='Incident flux (W / um2)')
    fig.suptitle('1 W per beam | same central incidence | separate absolute scales | measured demo is synthetic' if demo else
                 '1 W per beam | same central incidence | separate absolute scales')
    fig.savefig(output/'beam_definitions.png', dpi=150)
    plt.close(fig)
    result = dict(status='complete_with_physical_limits', cases=reports,
        scope='Incident field and propagating cutoff audit only, not a new FEM or scattering comparison. '
              'Intensity-only data use an explicit flat envelope phase; no phase retrieval is performed.')
    (output/'beam_definitions_report.json').write_text(json.dumps(result, indent=2), encoding='utf-8')
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', required=True)
    parser.add_argument('--output', required=True)
    parser.add_argument('--measurement')
    parser.add_argument('--phase-assumption', choices=('flat', 'reject'), default='flat')
    parser.add_argument('--native-check', action='store_true')
    args = parser.parse_args()
    print(json.dumps(run(args.source, args.output, args.measurement, args.phase_assumption, args.native_check), indent=2))


if __name__ == '__main__':
    main()
