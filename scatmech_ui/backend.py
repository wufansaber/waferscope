"""pySCATMECH model construction and far-field calculation."""

import numpy as np

try:
    from pySCATMECH.fresnel import Film, FilmStack, OpticalFunction
    from pySCATMECH.local import Local_BRDF_Model
    from pySCATMECH.mueller import Polarization

    SCATMECH_IMPORT_ERROR = None
except ModuleNotFoundError as error:
    if error.name not in ("pySCATMECH", "SCATPY"):
        raise
    Film = FilmStack = OpticalFunction = Local_BRDF_Model = Polarization = None
    SCATMECH_IMPORT_ERROR = error


def _build_particle_model(parameters):
    """Build a finite-film model where numerically meaningful."""
    wavelength_um = parameters.wavelength_nm / 1000.0
    wafer_index = complex(parameters.wafer_n, parameters.wafer_k)
    particle_index = complex(parameters.particle_n, parameters.particle_k)
    thickness = parameters.wafer_thickness_um
    optical_depth = 4.0 * np.pi * parameters.wafer_k * thickness / wavelength_um
    coherent_limit = 100.0 * wavelength_um

    model_parameters = {
        "lambda": wavelength_um,
        "radius": parameters.particle_diameter_nm / 2000.0,
        "sphere": particle_index,
        "type": 0,
    }

    if thickness <= coherent_limit and optical_depth < 60.0:
        wafer_stack = FilmStack(
            [Film(OpticalFunction(wafer_index), thickness=thickness)]
        )
        model_parameters["substrate"] = 1.0 + 0.0j
        model_parameters["stack"] = wafer_stack
        thickness_note = (
            "Finite coherent wafer layer: {:.6g} um on air backing."
        ).format(thickness)
    else:
        model_parameters["substrate"] = wafer_index
        thickness_note = (
            "Wafer is optically/geometrically thick at this wavelength; "
            "the backside is negligible and the wafer is modeled as semi-infinite."
        )

    model = Local_BRDF_Model("Bobbert_Vlieger_BRDF_Model", model_parameters)
    return model, thickness_note


def calculate_far_field(parameters, progress_callback=None):
    """Return intensity and full Stokes-derived polarization on a hemisphere."""
    if SCATMECH_IMPORT_ERROR is not None:
        raise RuntimeError(
            "pySCATMECH is unavailable in the active Python interpreter. "
            "Select a supported environment containing pySCATMECH and SCATPY."
        ) from SCATMECH_IMPORT_ERROR
    parameters.validate()
    model, thickness_note = _build_particle_model(parameters)
    step = parameters.angular_step_deg
    theta_deg = np.arange(0.0, 90.0, step)
    if theta_deg[-1] < 89.0:
        theta_deg = np.append(theta_deg, 89.0)
    phi_count = max(24, int(round(360.0 / step)))
    phi = np.linspace(0.0, 2.0 * np.pi, phi_count + 1)
    theta = np.deg2rad(theta_deg)
    phi_grid, theta_grid = np.meshgrid(phi, theta)

    shape = theta_grid.shape
    intensity = np.empty(shape, dtype=float)
    q = np.empty(shape, dtype=float)
    u = np.empty(shape, dtype=float)
    v = np.empty(shape, dtype=float)
    incident = Polarization(
        {"p": "p", "s": "s", "unpolarized": "u"}[parameters.polarization]
    )
    incidence = np.deg2rad(parameters.incidence_deg)

    for row in range(shape[0]):
        for column in range(shape[1]):
            mueller = model.MuellerDSC(
                incidence,
                theta_grid[row, column],
                phi_grid[row, column],
                coords="psps",
            )
            output = mueller @ incident
            intensity[row, column] = max(float(output[0]), 0.0)
            q[row, column] = float(output[1])
            u[row, column] = float(output[2])
            v[row, column] = float(output[3])
        if progress_callback is not None:
            progress_callback(int(round(100.0 * (row + 1) / shape[0])))

    if not np.isfinite(intensity).all() or np.nanmax(intensity) <= 0:
        raise RuntimeError("SCATMECH returned invalid far-field intensity values.")

    sin_theta = np.sin(theta_grid)
    cos_theta = np.cos(theta_grid)
    cos_phi = np.cos(phi_grid)
    sin_phi = np.sin(phi_grid)
    x = sin_theta * cos_phi
    y = sin_theta * sin_phi
    z = cos_theta

    safe_intensity = np.maximum(intensity, np.nanmax(intensity) * 1e-15)
    relative = safe_intensity / np.nanmax(safe_intensity)
    intensity_db = 10.0 * np.log10(relative)
    dolp = np.clip(np.sqrt(q**2 + u**2) / safe_intensity, 0.0, 1.0)
    docp = np.clip(v / safe_intensity, -1.0, 1.0)

    # Q is Is-Ip in SCATMECH's s-p basis. alpha=0 is s; alpha=pi/2 is p.
    alpha = 0.5 * np.arctan2(u, q)
    s_x = -sin_phi
    s_y = cos_phi
    s_z = np.zeros_like(s_x)
    p_x = cos_theta * cos_phi
    p_y = cos_theta * sin_phi
    p_z = -sin_theta
    e_x = np.cos(alpha) * s_x + np.sin(alpha) * p_x
    e_y = np.cos(alpha) * s_y + np.sin(alpha) * p_y
    e_z = np.cos(alpha) * s_z + np.sin(alpha) * p_z

    return {
        "parameters": parameters,
        "thickness_note": thickness_note,
        "theta": theta_grid,
        "phi": phi_grid,
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
    }
