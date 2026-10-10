"""Incident-only coherent beam, power normalization and rough-patch coverage."""
import json
from pathlib import Path

import numpy as np

from .finite_illumination import incident_fields, patch_coverage, power_area


def illumination_metadata(p, roughness=None):
    if p.beam_profile != 'gaussian':
        return dict(model='legacy_gaussian' if p.beam_waist_um else 'plane',
                    coherence='Single-frequency complex field; coherent amplitude sum',
                    finite_footprint=p.beam_waist_um > 0)
    from .beam_inputs import spectrum, truncation_report, wafer_dimensions
    directions, weights = spectrum(p)
    actual_length, actual_width = wafer_dimensions(p)
    report = dict(model='gaussian', coherence='Fully coherent angular-spectrum amplitude sum',
        finite_footprint=True, target_length_y_um=actual_length,
        target_width_x_um=actual_width, target_definition=p.beam_definition,
        beam_quadrature_order=p.beam_order, rays=len(directions),
        envelope_spectrum=truncation_report(p), normalization_area_um2=power_area(p),
        power_definition='Whole-wafer incident downward Poynting flux; continuous-spectrum Parseval '
                         'integral with '+('unit input-envelope peak E reference' if p.beam_definition in
                         ('wafer_rectangle', 'measured_wafer') else 'central reference flux normalization'),
        warnings=['Only propagating downward waves are included. Diffraction and grazing-angle '
                  'spectral clipping can deform the target envelope; inspect the actual beam map.'])
    if roughness is not None:
        coverage = patch_coverage(p, roughness.width_um)
        report.update(rough_patch_width_um=roughness.width_um,
                      incident_power_fraction_inside_rough_patch=coverage,
                      coverage_scope='Geometric patch before height taper; outside it the wafer is smooth')
        if coverage < .99:
            report['warnings'].append('The rough patch does not cover the complete illuminated wafer. '
                'Haze is from this modeled patch only; no area extrapolation is applied.')
    return report


def save_illumination(output, p, settings, roughness=None):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    report = illumination_metadata(p, roughness)
    length, width = settings.illumination_length_um or 2., settings.illumination_width_um or 2.
    if p.beam_profile == 'gaussian':
        from .beam_inputs import wafer_dimensions
        length, width = wafer_dimensions(p)
    x, y = np.linspace(-1.5*width, 1.5*width, 151), np.linspace(-1.5*length, 1.5*length, 151)
    xx, yy = np.meshgrid(x, y)
    points = np.stack((xx, yy, np.zeros_like(xx)), axis=-1)
    if p.beam_profile == 'gaussian':
        from .surface_beam import checked_surface_order, surface_fields_order
        order, checks, met = checked_surface_order(p, 3*max(length, width))
        _, e, h = surface_fields_order(x, y, p, order)
        report.update(display_field_quadrature_order=order, display_field_quadrature_checks=checks,
                      display_field_quadrature_status='checked' if met else 'target_not_met')
        if not met:
            report['warnings'].append('Displayed beam integral did not pass two stable refinements.')
    else:
        e, h = incident_fields(points, p)
    flux = -np.real(np.cross(e, h.conj()))[..., 2]*settings.irradiance_w_um2
    if p.beam_profile == 'gaussian':
        report['normalized_incident_power_w'] = settings.irradiance_w_um2*power_area(p)
    report['settings_power_w'] = settings.illumination_power_w
    np.savez_compressed(output/'illumination_fields.npz', x_um=x, y_um=y,
                        incident_e_unit=e, incident_hz0_unit=h, wafer_flux_w_um2=flux)
    fig, ax = plt.subplots(figsize=(7, 6), layout='constrained')
    image = ax.pcolormesh(x, y, flux, cmap='viridis', shading='auto', rasterized=True)
    if roughness is not None:
        from matplotlib.patches import Rectangle
        w = roughness.width_um
        ax.add_patch(Rectangle((-w/2, -w/2), w, w, fill=False, edgecolor='white', linewidth=1))
    ax.set(xlabel='Wafer X (um)', ylabel='Wafer Y (um)', aspect='equal',
           title='Incident-only coherent illumination\n{} | {} deg | {} nm'.format(
               report.get('target_definition', report['model']), p.incidence_deg, p.wavelength_nm))
    fig.colorbar(image, ax=ax, label='Downward incident flux (W / um2)')
    fig.savefig(output/'illumination.png', dpi=140)
    plt.close(fig)
    (output/'illumination_report.json').write_text(json.dumps(report, indent=2), encoding='utf-8')
    return report
