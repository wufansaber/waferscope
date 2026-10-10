"""Shared rough-surface geometry and absolute incident beam coverage for the UI."""
import json
from pathlib import Path

import numpy as np


def plot_surface(output, data, metadata, settings, audit):
    from .imaging_plotting import _pyplot
    from matplotlib.patches import Rectangle
    plt = _pyplot()
    output = Path(output)
    plan = metadata['plan']
    axis = data['whole_axis_um']
    local_width = metadata['roughness']['width_um']
    values = (data['whole_height_um']*data['whole_height_window']*1000,
              data['whole_incident_flux_unit']*settings.irradiance_w_um2)
    fig, axes = plt.subplots(1, 2, figsize=(12, 5), layout='constrained')
    for ax, value, title, unit, cmap in zip(axes, values,
            ('Shared whole rough surface', 'Incident-only coherent illumination'),
            ('Height (nm)', 'Downward flux (W / um2)'), ('coolwarm', 'viridis')):
        image = ax.pcolormesh(axis, axis, value, shading='auto', cmap=cmap)
        ax.add_patch(Rectangle((-local_width/2, -local_width/2), local_width, local_width,
                              fill=False, edgecolor='black', linewidth=1.2))
        ax.set(xlabel='Wafer X (um)', ylabel='Wafer Y (um)', aspect='equal', title=title)
        fig.colorbar(image, ax=ax, label=unit)
    fig.suptitle('Shared synthetic surface | whole width {:.4g} um | local FEM width {:.4g} um\n'
                 'Incident power {:.4g} W | height-window weighted coverage {:.3%}'.format(
                    plan['width_um'], local_width, settings.illumination_power_w,
                    plan['checked_taper_weighted_power_coverage']))
    fig.savefig(output/'illumination.png', dpi=140)
    plt.close(fig)
    step = plan['spacing_um']
    from .parameters import Parameters
    from .beam_inputs import wafer_dimensions, truncation_report
    p = Parameters(**metadata['parameters'])
    report = dict(model='hybrid_shared_surface_coherent_finite_beam',
        coherence='Fully coherent incoming angular-spectrum amplitudes; no random illumination phases',
        incident_total_power_w=settings.illumination_power_w,
        beam_definition=p.beam_definition, target_wafer_dimensions_um=list(wafer_dimensions(p)),
        envelope_spectrum=truncation_report(p),
        normalization_area_um2=settings.illumination_normalization_area_um2,
        whole_surface_plan=plan, local_fem_width_um=local_width,
        mother_surface=metadata['mother_surface'],
        geometric_input_power_coverage=float(np.sum(data['whole_incident_flux_unit'])*step**2/
                                              settings.illumination_normalization_area_um2),
        surface_convergence=audit,
        scope='Shared surface and beam for outer first-order SPM and the local FEM correction. '
              'Coverage is incident power, NOT collected scatter or an accuracy certificate.')
    (output/'illumination_report.json').write_text(json.dumps(report, indent=2), encoding='utf-8')
