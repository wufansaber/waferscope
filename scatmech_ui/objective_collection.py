"""Workflow for collecting far-field energy into ideal objective pupils."""

import numpy as np

from .objective_math import (
    dark_field_energy,
    format_tilt_list,
    local_pupil_directions,
    objective_far_field_view,
    pupil_polarization,
    reference_energy,
    relative_collected_power,
    sample_far_field,
    tilt_local_directions,
)


def calculate_objective_collection_example(parameters, source_far_field=None):
    """Map a far-field distribution into tilted ideal-objective pupils."""
    parameters.validate()
    reference_value = reference_energy(parameters, source_far_field)
    channel_pupils = [
        _objective_channel_pupil(
            parameters,
            tilt_deg,
            numerical_aperture,
            reference_value,
            source_far_field,
        )
        for tilt_deg, numerical_aperture in zip(
            parameters.tilt_degrees,
            parameters.channel_numerical_apertures,
        )
    ]
    far_field = objective_far_field_view(parameters, source_far_field)
    return _format_collection_result(parameters, channel_pupils, far_field)


def _objective_channel_pupil(
    parameters,
    tilt_deg,
    numerical_aperture,
    reference_value,
    source_far_field=None,
):
    local = local_pupil_directions(parameters, numerical_aperture)
    kx, ky, kz = tilt_local_directions(local, tilt_deg)
    aperture = local["aperture"] & (kz >= 0.0)
    if numerical_aperture == 0.0:
        aperture[:] = False
    phi = np.arctan2(ky, kx)
    far_field_per_sr, q, u = _sample_channel_field(source_far_field, kx, ky, kz)
    pupil_energy_density = far_field_per_sr / np.maximum(local["kz"], 1e-12)
    pupil_energy_density[~aperture] = np.nan
    intensity = np.clip(pupil_energy_density / reference_value, 0.0, 1.0)
    intensity_db = 10.0 * np.log10(
        np.maximum(np.nan_to_num(intensity, nan=0.0), 1e-15)
    )

    return {
        "tilt_deg": tilt_deg,
        "numerical_aperture": numerical_aperture,
        "pupil_radius_mm": parameters.focal_length_mm * numerical_aperture,
        "x": local["x"],
        "y": local["y"],
        "kx": kx,
        "ky": ky,
        "kz": kz,
        "local_kx": local["kx"],
        "local_ky": local["ky"],
        "theta": np.arccos(np.clip(kz, -1.0, 1.0)),
        "phi": phi,
        "intensity": intensity,
        "intensity_db": intensity_db,
        "pupil_energy_density": pupil_energy_density,
        "relative_collected_power": relative_collected_power(
            pupil_energy_density,
            reference_value,
            aperture,
        ),
        **pupil_polarization(
            intensity,
            far_field_per_sr,
            phi,
            kx,
            ky,
            kz,
            tilt_deg,
            aperture,
            q,
            u,
        ),
    }


def _sample_channel_field(source_far_field, kx, ky, kz):
    if source_far_field is None:
        return dark_field_energy(kx, ky), None, None
    sampled = sample_far_field(source_far_field, kx, ky, kz)
    return sampled["intensity"], sampled["q"], sampled["u"]


def _format_collection_result(parameters, channel_pupils, far_field):
    primary_pupil = channel_pupils[0]
    zero = np.zeros_like(primary_pupil["intensity"])
    return {
        **_primary_pupil_aliases(primary_pupil),
        **_far_field_aliases(far_field),
        "plot_kind": "objective_collection",
        "parameters": parameters,
        "thickness_note": _collection_note(parameters, far_field),
        "channel_pupils": channel_pupils,
        "q": zero,
        "u": zero,
        "v": zero,
        "dolp": zero,
        "docp": zero,
        "e_x": zero,
        "e_y": zero,
        "e_z": zero,
    }


def _primary_pupil_aliases(primary_pupil):
    return {
        "x": primary_pupil["x"],
        "y": primary_pupil["y"],
        "z": primary_pupil["intensity"],
        "kx": primary_pupil["kx"],
        "ky": primary_pupil["ky"],
        "kz": primary_pupil["kz"],
        "local_kx": primary_pupil["local_kx"],
        "local_ky": primary_pupil["local_ky"],
        "theta": primary_pupil["theta"],
        "phi": primary_pupil["phi"],
        "intensity": primary_pupil["intensity"],
        "intensity_db": primary_pupil["intensity_db"],
        "pupil_energy_density": primary_pupil["pupil_energy_density"],
        "relative_collected_power": primary_pupil["relative_collected_power"],
    }


def _far_field_aliases(far_field):
    aliases = {
        "collection_fraction": far_field["collection_fraction"],
        "energy_collection_fraction": far_field["energy_collection_fraction"],
        "geometric_collection_fraction": far_field[
            "geometric_collection_fraction"
        ],
        "collected_energy": far_field["collected_energy"],
        "hemisphere_energy": far_field["hemisphere_energy"],
        "geometric_collected_volume": far_field["geometric_collected_volume"],
        "hemisphere_geometric_volume": far_field[
            "hemisphere_geometric_volume"
        ],
        "far_x": far_field["x"],
        "far_y": far_field["y"],
        "far_z": far_field["z"],
        "far_intensity": far_field["intensity"],
        "far_intensity_db": far_field["intensity_db"],
        "far_collection_mask": far_field["collection_mask"],
        "far_channel_masks": far_field["channel_masks"],
        "channel_fractions": far_field["channel_fractions"],
        "energy_channel_fractions": far_field["energy_channel_fractions"],
        "geometric_channel_fractions": far_field[
            "geometric_channel_fractions"
        ],
        "channel_collected_energies": far_field["channel_collected_energies"],
        "channel_geometric_volumes": far_field["channel_geometric_volumes"],
        "objective_axis": far_field["objective_axis"],
        "objective_axes": far_field["objective_axes"],
        "cone_boundary": far_field["cone_boundary"],
        "cone_boundaries": far_field["cone_boundaries"],
    }
    for key in (
        "incidence_plane",
        "illumination_rotation_deg",
        "scatmech_azimuth_offset_deg",
        "incident_direction",
        "incident_source_direction",
        "incident_projection_direction",
    ):
        if key in far_field:
            aliases[key] = far_field[key]
    return aliases


def _collection_note(parameters, far_field):
    source_note = (
        "SCATMECH far-field"
        if far_field["source"] == "scatmech"
        else "synthetic dark-field"
    )
    na_text = "/".join(
        "{:.3g}".format(value)
        for value in parameters.channel_numerical_apertures
    )
    radius_text = "/".join(
        "{:.3g}".format(parameters.focal_length_mm * value)
        for value in parameters.channel_numerical_apertures
    )
    collection_note = (
        "Ideal single-lens objective collection from {}: channel NA={}, "
        "f={:.3g} mm, tilts={}, exit-pupil radii={} mm, "
        "union geometric volume={:.2f}%, union energy={:.2f}%."
    ).format(
        source_note,
        na_text,
        parameters.focal_length_mm,
        format_tilt_list(parameters.tilt_degrees),
        radius_text,
        100.0 * far_field["geometric_collection_fraction"],
        100.0 * far_field["collection_fraction"],
    )
    source_notes = [
        far_field.get("thickness_note", ""),
        far_field.get("model_note", ""),
    ]
    return " ".join(
        [collection_note] + [note for note in source_notes if note]
    )
