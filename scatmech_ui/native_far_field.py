"""Native pySCATMECH far-field model construction and sampling."""

import numpy as np

from .coordinate_geometry import (
    INCIDENCE_PLANE,
    ILLUMINATION_ROTATION_DEG,
    SCATMECH_AZIMUTH_OFFSET_DEG,
    incident_projection_direction,
    incident_propagation_direction,
    incident_source_direction,
    scatmech_to_global_xy,
)
from .surface_model import (
    Polarization,
    SCATMECH_IMPORT_ERROR,
    build_scattering_model,
)


def calculate_far_field(parameters, progress_callback=None):
    """Return intensity and Stokes-derived polarization on a hemisphere."""
    if SCATMECH_IMPORT_ERROR is not None:
        raise RuntimeError(
            "pySCATMECH is unavailable in the active Python interpreter. "
            "Select a supported environment containing pySCATMECH and SCATPY."
        ) from SCATMECH_IMPORT_ERROR

    parameters.validate()
    model_bundle = build_scattering_model(parameters)
    theta_grid, phi_grid = _angular_grid(parameters.angular_step_deg)
    intensity, q, u, v = _evaluate_mueller_grid(
        model_bundle,
        parameters,
        theta_grid,
        phi_grid,
        progress_callback,
    )
    return _format_far_field_result(
        parameters,
        model_bundle.thickness_note,
        theta_grid,
        phi_grid,
        intensity,
        q,
        u,
        v,
        model_note=model_bundle.model_note,
    )


def _angular_grid(step_deg):
    theta_deg = np.arange(0.0, 90.0, step_deg)
    if theta_deg[-1] < 89.0:
        theta_deg = np.append(theta_deg, 89.0)
    phi_count = max(24, int(round(360.0 / step_deg)))
    phi = np.linspace(0.0, 2.0 * np.pi, phi_count + 1)
    theta = np.deg2rad(theta_deg)
    phi_grid, theta_grid = np.meshgrid(phi, theta)
    _validate_angular_grid(theta_grid, phi_grid)
    return theta_grid, phi_grid


def _validate_angular_grid(theta_grid, phi_grid):
    if np.nanmax(theta_grid) > 0.5 * np.pi or np.nanmin(theta_grid) < 0.0:
        raise ValueError("theta grid must be within the upper hemisphere [0, pi/2].")
    if np.nanmax(phi_grid) > 2.0 * np.pi or np.nanmin(phi_grid) < 0.0:
        raise ValueError("phi grid must be within [0, 2*pi].")


def _evaluate_mueller_grid(
    model_bundle,
    parameters,
    theta_grid,
    phi_grid,
    progress_callback,
):
    shape = theta_grid.shape
    intensity = np.empty(shape, dtype=float)
    q = np.empty(shape, dtype=float)
    u = np.empty(shape, dtype=float)
    v = np.empty(shape, dtype=float)
    incident = Polarization(
        {"p": "p", "s": "s", "unpolarized": "u"}[parameters.polarization]
    )
    incidence = np.deg2rad(parameters.incidence_deg)
    incidence_cosine = abs(np.cos(incidence))

    for row in range(shape[0]):
        for column in range(shape[1]):
            scattering_theta = theta_grid[row, column]
            if model_bundle.output_kind == "dsc":
                mueller = model_bundle.model.MuellerDSC(
                    incidence,
                    scattering_theta,
                    phi_grid[row, column],
                    coords="psps",
                )
                output_scale = 1.0
            else:
                mueller = model_bundle.model.MuellerBRDF(
                    incidence,
                    scattering_theta,
                    phi_grid[row, column],
                    coords="psps",
                )
                output_scale = (
                    model_bundle.effective_area_um2
                    * incidence_cosine
                    * abs(np.cos(scattering_theta))
                )
            output = mueller @ incident
            intensity[row, column] = max(float(output[0]) * output_scale, 0.0)
            q[row, column] = float(output[1]) * output_scale
            u[row, column] = float(output[2]) * output_scale
            v[row, column] = float(output[3]) * output_scale
        if progress_callback is not None:
            progress_callback(int(round(100.0 * (row + 1) / shape[0])))

    if not np.isfinite(intensity).all() or np.nanmax(intensity) <= 0:
        raise RuntimeError("SCATMECH returned invalid far-field intensity values.")
    return intensity, q, u, v


def _format_far_field_result(
    parameters,
    thickness_note,
    theta_grid,
    phi_grid,
    intensity,
    q,
    u,
    v,
    model_note="",
):
    sin_theta = np.sin(theta_grid)
    cos_theta = np.cos(theta_grid)
    cos_phi = np.cos(phi_grid)
    sin_phi = np.sin(phi_grid)
    local_x = sin_theta * cos_phi
    local_y = sin_theta * sin_phi
    x, y = scatmech_to_global_xy(local_x, local_y)
    z = cos_theta

    safe_intensity = np.maximum(intensity, np.nanmax(intensity) * 1e-15)
    relative = safe_intensity / np.nanmax(safe_intensity)
    intensity_db = 10.0 * np.log10(relative)
    dolp = np.clip(np.sqrt(q**2 + u**2) / safe_intensity, 0.0, 1.0)
    docp = np.clip(v / safe_intensity, -1.0, 1.0)

    alpha = 0.5 * np.arctan2(u, q)
    local_s_x = -sin_phi
    local_s_y = cos_phi
    s_x, s_y = scatmech_to_global_xy(local_s_x, local_s_y)
    s_z = np.zeros_like(s_x)
    local_p_x = cos_theta * cos_phi
    local_p_y = cos_theta * sin_phi
    p_x, p_y = scatmech_to_global_xy(local_p_x, local_p_y)
    p_z = -sin_theta
    e_x = np.cos(alpha) * s_x + np.sin(alpha) * p_x
    e_y = np.cos(alpha) * s_y + np.sin(alpha) * p_y
    e_z = np.cos(alpha) * s_z + np.sin(alpha) * p_z
    incident_direction = incident_propagation_direction(parameters.incidence_deg)

    return {
        "parameters": parameters,
        "thickness_note": thickness_note,
        "model_note": model_note,
        "theta": theta_grid,
        "phi": phi_grid,
        "scatmech_phi": phi_grid,
        "global_phi": np.mod(
            phi_grid + np.deg2rad(SCATMECH_AZIMUTH_OFFSET_DEG),
            2.0 * np.pi,
        ),
        "x": x,
        "y": y,
        "z": z,
        "intensity": intensity,
        "intensity_db": intensity_db,
        "q": q,
        "u": u,
        "v": v,
        "dolp": dolp,
        "docp": docp,
        "e_x": e_x,
        "e_y": e_y,
        "e_z": e_z,
        "incidence_plane": INCIDENCE_PLANE,
        "illumination_rotation_deg": ILLUMINATION_ROTATION_DEG,
        "scatmech_azimuth_offset_deg": SCATMECH_AZIMUTH_OFFSET_DEG,
        "incident_direction": np.asarray(incident_direction, dtype=float),
        "incident_source_direction": np.asarray(
            incident_source_direction(parameters.incidence_deg),
            dtype=float,
        ),
        "incident_projection_direction": np.asarray(
            incident_projection_direction(parameters.incidence_deg),
            dtype=float,
        ),
    }
