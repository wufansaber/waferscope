"""Translate BV sphere-centered amplitudes to the top-wafer phase origin."""
import numpy as np


def sphere_to_wafer_phase(directions, wavelength_um, center_height_um, incidence_deg):
    directions = np.asarray(directions, dtype=float)
    if directions.shape[-1:] != (3,) or not np.isfinite(directions).all():
        raise ValueError('Finite three-component observation directions required')
    if (not np.isfinite([wavelength_um, center_height_um, incidence_deg]).all()
            or wavelength_um <= 0 or center_height_um < 0 or not 0 <= incidence_deg < 90):
        raise ValueError('Invalid sphere-to-wafer phase-reference parameters')
    # Incident k_z=-cos(theta_i); outgoing translation contributes -k_s dot r_c.
    phase = -2*np.pi/wavelength_um*center_height_um*(
        directions[..., 2]+np.cos(np.deg2rad(incidence_deg)))
    return np.exp(1j*phase)
