"""Shared ideal aplanatic vector propagation, independent of UI and scattering backend."""
import numpy as np


def pupil_geometry(na, tilt_deg, samples):
    # An even midpoint grid avoids the s/p singularity at the optical axis.
    du = 2.0 * (na if na > 0 else 0.1) / samples
    axis = (np.arange(samples) - (samples - 1) / 2.0) * du
    u, v = np.meshgrid(axis, axis)
    w = np.sqrt(np.maximum(1.0 - u*u - v*v, 1e-12))
    local = np.stack((u, v, w), axis=-1)
    t = np.deg2rad(tilt_deg)
    rotation = np.array(((np.cos(t), 0, np.sin(t)),
                         (0, 1, 0), (-np.sin(t), 0, np.cos(t))))
    global_k = local @ rotation.T
    mask = (u*u + v*v < na*na) & (global_k[..., 2] > 1e-8)
    return local, global_k, rotation, mask, du


def collimated_pupil(jones, global_k, local_k, rotation):
    rho = np.hypot(global_k[..., 0], global_k[..., 1])
    s = np.stack((-global_k[..., 1], global_k[..., 0],
                  np.zeros_like(rho)), axis=-1) / np.maximum(rho[..., None], 1e-15)
    # SCATMECH uses local s=(0,1,0) on axis, rotated to global (-1,0,0).
    s[rho < 1e-12] = (-1, 0, 0)
    p = np.cross(global_k, s)
    basis = np.stack((s, p), axis=-1)
    field = np.einsum("...ci,...ij->...cj", basis, jones)
    field = np.einsum("dc,...ci->...di", rotation.T, field)
    # Minimal ray rotation k -> z, preserving both polarizations and power.
    xy = field[..., :2, :] - (
        local_k[..., :2, None] * field[..., 2:3, :]
        / (1.0 + local_k[..., 2, None, None])
    )
    return xy / np.sqrt(local_k[..., 2, None, None])


def propagate_pupil(pupil, local_k, du, wavelength_um, parameters):
    """Return vector image modes and H*Z0 modes with exact axial flux weights."""
    m = parameters.magnification
    k_image = np.empty_like(local_k)
    k_image[..., :2] = -local_k[..., :2] / m
    propagating = np.sum(k_image[..., :2]**2, axis=-1) < 1
    if np.any(~propagating & np.any(abs(pupil) > 0, axis=(-2, -1))):
        raise ValueError('Nonzero pupil field contains non-propagating image directions')
    # Empty square-grid corners can lie outside the unit disk at high image NA.
    # Give those zero-field cells a finite dummy direction before vector rotation.
    k_image[..., :2] = np.where(propagating[..., None], k_image[..., :2], 0)
    k_image[..., 2] = np.sqrt(1 - np.sum(k_image[..., :2]**2, axis=-1))
    uv = k_image[..., :2]
    dot = np.einsum("...c,...ci->...i", uv, pupil)
    e = np.zeros(pupil.shape[:-2] + (3, 2), dtype=complex)
    e[..., :2, :] = pupil - uv[..., :, None] * dot[..., None, :] / (
        1 + k_image[..., 2, None, None]
    )
    e[..., 2, :] = -dot
    e /= np.sqrt(k_image[..., 2, None, None])
    phase = np.exp(2j*np.pi/wavelength_um * parameters.image_defocus_um
                   * (k_image[..., 2] - 1))
    e *= phase[..., None, None]
    h = np.moveaxis(np.cross(k_image[..., None, :], np.moveaxis(e, -2, -1)), -1, -2)
    n = pupil.shape[0]
    size = n * parameters.padding
    pad = (size - n) // 2

    def transform(values):
        padded = np.pad(values, ((pad, pad), (pad, pad), (0, 0), (0, 0)))
        transformed = np.fft.fftshift(np.fft.fft2(
            np.fft.ifftshift(padded, axes=(0, 1)), axes=(0, 1)), axes=(0, 1))
        return transformed * du**2 / (wavelength_um * m)

    axis = np.fft.fftshift(np.fft.fftfreq(size, d=du)) * wavelength_um * m
    # Midpoint input grid is offset +du/2 from the FFT origin.
    frequency = axis / (wavelength_um * m)
    correction = np.exp(-1j*np.pi*du*(frequency[:, None] + frequency[None, :]))
    return axis, transform(e)*correction[..., None, None], transform(h)*correction[..., None, None]
