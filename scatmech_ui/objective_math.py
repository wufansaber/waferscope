"""Math helpers for ideal-objective far-field collection."""

import numpy as np


def reference_energy(parameters, source_far_field):
    if source_far_field is not None:
        return max(float(np.nanmax(source_far_field["intensity"])), 1e-15)
    maximum_na = max(parameters.channel_numerical_apertures)
    return 4.0 / max(np.sqrt(1.0 - maximum_na**2), 1e-6)


def local_pupil_directions(parameters, numerical_aperture):
    pupil_radius_mm = parameters.focal_length_mm * numerical_aperture
    coordinate = np.linspace(-pupil_radius_mm, pupil_radius_mm, parameters.grid_size)
    x, y = np.meshgrid(coordinate, coordinate)
    sin_theta = np.clip(np.hypot(x, y) / parameters.focal_length_mm, 0.0, 1.0)
    aperture = sin_theta <= numerical_aperture
    return {
        "x": x,
        "y": y,
        "kx": x / parameters.focal_length_mm,
        "ky": y / parameters.focal_length_mm,
        "kz": np.sqrt(np.maximum(1.0 - sin_theta**2, 1e-12)),
        "aperture": aperture,
    }


def tilt_local_directions(local, tilt_deg):
    tilt = np.deg2rad(tilt_deg)
    kx = np.cos(tilt) * local["kx"] + np.sin(tilt) * local["kz"]
    ky = local["ky"]
    kz = -np.sin(tilt) * local["kx"] + np.cos(tilt) * local["kz"]
    return kx, ky, kz


def relative_collected_power(pupil_energy_density, reference_value, aperture):
    collected_power = float(np.nansum(pupil_energy_density))
    reference_power = float(reference_value * np.count_nonzero(aperture))
    return collected_power / max(reference_power, 1e-12)


def pupil_polarization(intensity, far_field_per_sr, phi, kx, ky, aperture, q, u):
    if q is None or u is None:
        high_angle = np.clip(np.hypot(kx, ky), 0.0, 1.0)
        dolp = np.clip(0.20 + 0.70 * high_angle**2, 0.0, 0.95)
        aolp = phi + 0.5 * np.pi + 0.18 * np.sin(2.0 * phi)
        q = np.nan_to_num(intensity, nan=0.0) * dolp * np.cos(2.0 * aolp)
        u = np.nan_to_num(intensity, nan=0.0) * dolp * np.sin(2.0 * aolp)
    else:
        safe_far_field = np.maximum(far_field_per_sr, 1e-15)
        dolp = np.clip(np.sqrt(q**2 + u**2) / safe_far_field, 0.0, 1.0)
        aolp = 0.5 * np.arctan2(u, q)

    dolp[~aperture] = np.nan
    return {
        "dolp": dolp,
        "aolp": aolp,
        "pol_x": np.cos(aolp),
        "pol_y": np.sin(aolp),
        "q": q,
        "u": u,
    }


def objective_far_field_view(parameters, source_far_field=None):
    far_field = (
        synthetic_far_field(parameters)
        if source_far_field is None
        else scatmech_far_field(source_far_field)
    )
    return {
        **far_field,
        **collection_masks_and_energy(parameters, far_field),
    }


def synthetic_far_field(parameters):
    theta = np.linspace(0.0, 0.5 * np.pi, parameters.grid_size)
    phi = np.linspace(0.0, 2.0 * np.pi, parameters.grid_size)
    phi_grid, theta_grid = np.meshgrid(phi, theta)
    sin_theta = np.sin(theta_grid)
    x = sin_theta * np.cos(phi_grid)
    y = sin_theta * np.sin(phi_grid)
    z = np.cos(theta_grid)
    energy = dark_field_energy(x, y)
    maximum_na = max(parameters.channel_numerical_apertures)
    reference = 4.0 / max(np.sqrt(1.0 - maximum_na**2), 1e-6)
    intensity = np.clip(energy / reference, 0.0, 1.0)
    return far_field_view_dict("synthetic", theta_grid, x, y, z, energy, intensity)


def scatmech_far_field(source_far_field):
    energy = source_far_field["intensity"]
    reference = max(float(np.nanmax(energy)), 1e-15)
    intensity = np.clip(energy / reference, 0.0, 1.0)
    far_field = far_field_view_dict(
        "scatmech",
        source_far_field["theta"],
        source_far_field["x"],
        source_far_field["y"],
        source_far_field["z"],
        energy,
        intensity,
    )
    for key in (
        "incidence_plane",
        "illumination_rotation_deg",
        "incident_direction",
        "incident_source_direction",
        "incident_projection_direction",
    ):
        if key in source_far_field:
            far_field[key] = source_far_field[key]
    return far_field


def far_field_view_dict(source, theta, x, y, z, energy, intensity):
    return {
        "source": source,
        "theta": theta,
        "x": x,
        "y": y,
        "z": z,
        "energy": energy,
        "intensity": intensity,
        "intensity_db": 10.0 * np.log10(np.maximum(intensity, 1e-15)),
    }


def collection_masks_and_energy(parameters, far_field):
    mask_geometry = objective_collection_masks(
        parameters,
        far_field["x"],
        far_field["y"],
        far_field["z"],
    )
    energy_result = collection_energy(
        far_field["theta"],
        far_field["energy"],
        mask_geometry["collection_mask"],
        mask_geometry["channel_masks"],
    )
    return {
        **mask_geometry,
        **energy_result,
        **collection_geometry(parameters),
    }


def objective_collection_masks(parameters, x, y, z):
    objective_axes = []
    cone_boundaries = []
    channel_masks = []
    collection_mask = np.zeros_like(x, dtype=bool)
    for tilt_deg, numerical_aperture in zip(
        parameters.tilt_degrees,
        parameters.channel_numerical_apertures,
    ):
        axis, boundary = objective_axis_and_boundary(
            tilt_deg,
            numerical_aperture,
        )
        objective_axes.append(axis)
        if numerical_aperture == 0.0:
            cone_boundaries.append(None)
            channel_masks.append(np.zeros_like(collection_mask))
            continue

        cone_half_angle = np.arcsin(numerical_aperture)
        cone_boundaries.append(boundary)
        direction_dot_axis = (
            x * axis[0]
            + y * axis[1]
            + z * axis[2]
        )
        channel_mask = direction_dot_axis >= np.cos(cone_half_angle)
        channel_masks.append(channel_mask)
        collection_mask |= channel_mask

    return {
        "collection_mask": collection_mask,
        "channel_masks": channel_masks,
        "objective_axis": objective_axes[0],
        "objective_axes": objective_axes,
        "cone_boundary": cone_boundaries[0],
        "cone_boundaries": cone_boundaries,
    }


def collection_energy(
    theta_grid,
    energy,
    collection_mask,
    channel_masks,
):
    solid_angle_weight = np.sin(theta_grid)
    hemisphere_energy = float(np.sum(energy * solid_angle_weight))
    channel_collected_energies = [
        float(np.sum(energy[mask] * solid_angle_weight[mask]))
        for mask in channel_masks
    ]
    collected_energy = float(
        np.sum(energy[collection_mask] * solid_angle_weight[collection_mask])
    )
    return {
        "collection_mask": collection_mask,
        "channel_masks": channel_masks,
        "channel_collected_energies": channel_collected_energies,
        "channel_fractions": [
            energy_value / max(hemisphere_energy, 1e-12)
            for energy_value in channel_collected_energies
        ],
        "collection_fraction": collected_energy / max(hemisphere_energy, 1e-12),
        "energy_channel_fractions": [
            energy_value / max(hemisphere_energy, 1e-12)
            for energy_value in channel_collected_energies
        ],
        "energy_collection_fraction": collected_energy
        / max(hemisphere_energy, 1e-12),
        "collected_energy": collected_energy,
        "hemisphere_energy": hemisphere_energy,
    }


def collection_geometry(parameters):
    """Calculate unit-hemisphere cone volumes without intensity weighting."""
    theta_count = 900
    phi_count = 3600
    theta_step = 0.5 * np.pi / theta_count
    phi_step = 2.0 * np.pi / phi_count
    theta = (np.arange(theta_count) + 0.5)[:, None] * theta_step
    phi = (np.arange(phi_count) + 0.5)[None, :] * phi_step
    sin_theta = np.sin(theta)
    x = sin_theta * np.cos(phi)
    y = sin_theta * np.sin(phi)
    z = np.cos(theta)
    geometry = objective_collection_masks(parameters, x, y, z)
    solid_angle_weight = np.broadcast_to(sin_theta, x.shape)
    solid_angle_scale = theta_step * phi_step
    hemisphere_solid_angle = float(np.sum(solid_angle_weight) * solid_angle_scale)
    channel_solid_angles = [
        float(np.sum(solid_angle_weight[mask]) * solid_angle_scale)
        for mask in geometry["channel_masks"]
    ]
    collected_solid_angle = float(
        np.sum(solid_angle_weight[geometry["collection_mask"]])
        * solid_angle_scale
    )
    return {
        "geometric_channel_fractions": [
            solid_angle / hemisphere_solid_angle
            for solid_angle in channel_solid_angles
        ],
        "geometric_collection_fraction": collected_solid_angle
        / hemisphere_solid_angle,
        "channel_geometric_volumes": [
            solid_angle / 3.0 for solid_angle in channel_solid_angles
        ],
        "geometric_collected_volume": collected_solid_angle / 3.0,
        "hemisphere_geometric_volume": hemisphere_solid_angle / 3.0,
    }


def objective_axis_and_boundary(tilt_deg, numerical_aperture):
    tilt = np.deg2rad(tilt_deg)
    axis = np.array([np.sin(tilt), 0.0, np.cos(tilt)])
    azimuth = np.linspace(0.0, 2.0 * np.pi, 241)
    local_x = numerical_aperture * np.cos(azimuth)
    local_y = numerical_aperture * np.sin(azimuth)
    local_z = np.full_like(
        local_x,
        np.sqrt(max(1.0 - numerical_aperture**2, 0.0)),
    )
    boundary_x = np.cos(tilt) * local_x + np.sin(tilt) * local_z
    boundary_y = local_y
    boundary_z = -np.sin(tilt) * local_x + np.cos(tilt) * local_z
    return axis, (boundary_x, boundary_y, boundary_z)


def sample_far_field(source_far_field, kx, ky, kz):
    theta = np.arccos(np.clip(kz, -1.0, 1.0))
    phi = np.mod(np.arctan2(ky, kx), 2.0 * np.pi)
    theta_axis = source_far_field["theta"][:, 0]
    phi_axis = source_far_field["phi"][0, :]
    theta_flat = theta.ravel()
    phi_flat = phi.ravel()

    theta_index = np.searchsorted(theta_axis, theta_flat, side="right") - 1
    phi_index = np.searchsorted(phi_axis, phi_flat, side="right") - 1
    theta_index = np.clip(theta_index, 0, len(theta_axis) - 2)
    phi_index = np.clip(phi_index, 0, len(phi_axis) - 2)

    theta_weight = _axis_weight(theta_flat, theta_axis, theta_index)
    phi_weight = _axis_weight(phi_flat, phi_axis, phi_index)

    return {
        "intensity": _interpolate_field(
            source_far_field["intensity"],
            theta_index,
            phi_index,
            theta_weight,
            phi_weight,
            theta.shape,
        ),
        "q": _interpolate_field(
            source_far_field["q"],
            theta_index,
            phi_index,
            theta_weight,
            phi_weight,
            theta.shape,
        ),
        "u": _interpolate_field(
            source_far_field["u"],
            theta_index,
            phi_index,
            theta_weight,
            phi_weight,
            theta.shape,
        ),
    }


def _axis_weight(values, axis_values, index):
    axis0 = axis_values[index]
    axis1 = axis_values[index + 1]
    return (values - axis0) / np.maximum(axis1 - axis0, 1e-12)


def _interpolate_field(values, theta_index, phi_index, theta_weight, phi_weight, shape):
    v00 = values[theta_index, phi_index]
    v01 = values[theta_index, phi_index + 1]
    v10 = values[theta_index + 1, phi_index]
    v11 = values[theta_index + 1, phi_index + 1]
    interpolated = (
        (1.0 - theta_weight) * (1.0 - phi_weight) * v00
        + (1.0 - theta_weight) * phi_weight * v01
        + theta_weight * (1.0 - phi_weight) * v10
        + theta_weight * phi_weight * v11
    )
    return interpolated.reshape(shape)


def dark_field_energy(kx, ky):
    high_angle = np.clip(np.hypot(kx, ky), 0.0, 1.0)
    phi = np.arctan2(ky, kx)
    dark_field_growth = 0.02 + 0.18 * high_angle**2 + 2.6 * high_angle**6
    forward_lobe = 0.55 * np.exp(
        -(((kx - 0.72) / 0.20) ** 2 + (ky / 0.28) ** 2)
    )
    weak_side_lobe = 0.18 * np.exp(
        -(((kx + 0.45) / 0.24) ** 2 + ((ky - 0.34) / 0.22) ** 2)
    )
    azimuthal_texture = 0.15 * (1.0 + np.cos(4.0 * phi)) * high_angle**2
    return (dark_field_growth + forward_lobe + weak_side_lobe) * (
        1.0 + azimuthal_texture
    )


def format_tilt_list(tilt_degrees):
    return "/".join("{:.3g} deg".format(tilt_deg) for tilt_deg in tilt_degrees)
