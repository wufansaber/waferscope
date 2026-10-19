"""Coherent particle imaging through independent ideal aplanatic objectives.

JonesDSC uses SCATMECH's (s, p) basis, with p = k cross s. Fields are
flux-normalized, not SI electric fields. All integration lengths are in um.
"""

from dataclasses import asdict, dataclass, field
import json
from pathlib import Path

import numpy as np
from optical_propagation import pupil_geometry, collimated_pupil, propagate_pupil

from .coordinate_geometry import incidence_metadata
from .surface_model import build_scattering_model
from .detector import DetectorParameters, calculate_detector


@dataclass(frozen=True)
class ImagingParameters:
    magnification: float = 20.0
    irradiance_mw_um2: float = 0.001
    image_defocus_um: float = 0.0
    transmission: float = 1.0
    analyzer: str = "none"
    pupil_samples: int = 64
    padding: int = 4
    multipole_order: int = 0
    detector: DetectorParameters = field(default_factory=DetectorParameters)

    def validate(self):
        self.detector.validate()
        values = (self.magnification, self.irradiance_mw_um2,
                  self.image_defocus_um, self.transmission)
        if not np.isfinite(values).all():
            raise ValueError("Imaging parameters must be finite.")
        if self.magnification < 1 or self.irradiance_mw_um2 <= 0:
            raise ValueError("Magnification must be >= 1 and irradiance positive.")
        if not 0 <= self.transmission <= 1:
            raise ValueError("Power transmission must be between 0 and 1.")
        if self.analyzer not in ("none", "x", "y", "45", "135"):
            raise ValueError("Unknown pupil analyzer.")
        if self.pupil_samples not in (32, 64, 96, 128):
            raise ValueError("Pupil samples must be 32, 64, 96, or 128.")
        if self.padding not in (2, 4, 8):
            raise ValueError("FFT padding must be 2, 4, or 8.")
        if not isinstance(self.multipole_order, int) or not 0 <= self.multipole_order <= 60:
            raise ValueError("Multipole order must be an integer from 0 (auto) to 60.")




def sample_jones_pupil(model, incidence, global_k, mask, progress=None):
    jones = np.zeros(mask.shape + (2, 2), dtype=complex)
    theta = np.arccos(np.clip(global_k[..., 2], -1, 1))
    # Global axes are SCATMECH axes rotated by +90 degrees around Z.
    phi = np.arctan2(-global_k[..., 0], global_k[..., 1])
    for row in range(mask.shape[0]):
        for col in np.flatnonzero(mask[row]):
            jones[row, col] = model.JonesDSC(
                incidence, theta[row, col], phi[row, col], coords="psps"
            )
        if progress:
            progress((row + 1) / mask.shape[0])
    if not np.isfinite(jones).all():
        raise RuntimeError("Non-finite SCATMECH Jones pupil.")
    return jones




def mode_weights(polarization):
    return {"s": np.array((1., 0.)), "p": np.array((0., 1.)),
            "unpolarized": np.array((0.5, 0.5))}[polarization]


def _channel(model, simulation, objective, imaging, na, tilt, progress):
    local, global_k, rotation, mask, du = pupil_geometry(na, tilt, imaging.pupil_samples)
    jones = sample_jones_pupil(model, np.deg2rad(simulation.incidence_deg),
                               global_k, mask, progress)
    pupil = collimated_pupil(jones, global_k, local, rotation)
    weights = mode_weights(simulation.polarization)
    collected = float(np.sum(np.abs(pupil)**2 * weights) * du**2)
    if imaging.analyzer != "none":
        angle = {"x": 0, "y": 90, "45": 45, "135": 135}[imaging.analyzer]
        a = np.array((np.cos(np.deg2rad(angle)), np.sin(np.deg2rad(angle))))
        projection = np.einsum("c,...ci->...i", a, pupil)
        pupil = a[None, None, :, None] * projection[..., None, :]
    pupil *= np.sqrt(imaging.transmission)
    transmitted = float(np.sum(np.abs(pupil)**2 * weights) * du**2)
    axis, e, h = propagate_pupil(pupil, local, du, simulation.wavelength_nm/1000, imaging)
    flux_modes = np.real(e[..., 0, :] * h[..., 1, :].conj()
                         - e[..., 1, :] * h[..., 0, :].conj())
    flux = np.sum(flux_modes * weights, axis=-1)
    image_dsc = float(flux.sum() * (axis[1]-axis[0])**2)
    error = abs(image_dsc-transmitted)/transmitted if transmitted > 0 else 0.0
    irradiance = imaging.irradiance_mw_um2 * 1e-3
    edge = np.ones(flux.shape, dtype=bool)
    border = max(1, len(axis)//10)
    edge[border:-border, border:-border] = False
    edge_fraction = float(np.abs(flux[edge]).sum()/max(np.abs(flux).sum(), 1e-300))
    specular = np.array((0, np.sin(np.deg2rad(simulation.incidence_deg)),
                         np.cos(np.deg2rad(simulation.incidence_deg))))
    admits_specular = na > 0 and float(specular @ rotation[:, 2]) >= np.sqrt(1-na*na)
    return dict(na=na, tilt_deg=tilt, image_axis_um=axis,
                pupil_axis_mm=local[0, :, 0]*objective.focal_length_mm,
                pupil_mask=mask, jones_dsc_um=jones, pupil_modes=pupil,
                image_e_modes=e, image_hz0_modes=h,
                irradiance_w_um2=flux*irradiance,
                collected_dsc_um2=collected, transmitted_dsc_um2=transmitted,
                image_integral_dsc_um2=image_dsc, power_w=image_dsc*irradiance,
                flux_relative_error=error, edge_fraction=edge_fraction,
                admits_specular=bool(admits_specular))


def calculate_imaging(simulation, objective, imaging, progress_callback=None):
    simulation.validate()
    objective.validate()
    imaging.validate()
    if simulation.roughness_enabled:
        raise ValueError("Coherent imaging requires roughness disabled: a mean BRDF "
                         "does not specify spatial phase correlations.")
    bundle = build_scattering_model(simulation)
    bundle.model.setParameters(lmax=imaging.multipole_order, order=-1, Norm_Inc_Approx=0)
    channels = []
    streams = np.random.SeedSequence(imaging.detector.seed).spawn(3)
    for index, (na, tilt) in enumerate(zip(objective.channel_numerical_apertures,
                                         objective.tilt_degrees)):
        def progress(fraction):
            if progress_callback:
                progress_callback(int(100*(index+fraction)/3))
        channel = _channel(bundle.model, simulation, objective, imaging, na, tilt, progress)
        if imaging.detector.enabled:
            channel["detector"] = calculate_detector(
                channel, simulation.wavelength_nm, imaging.detector,
                np.random.default_rng(streams[index]), sensor_index=index)
        channels.append(channel)
    return dict(plot_kind="coherent_imaging", parameters=simulation,
                objective_parameters=objective, imaging_parameters=imaging,
                tube_focal_length_mm=objective.focal_length_mm*imaging.magnification,
                channels=channels, model_note=bundle.model_note,
                thickness_note=bundle.thickness_note,
                near_field_available=False, **incidence_metadata(simulation.incidence_deg))


def export_imaging(result, directory, figure, detector_figure=None):
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    figure.savefig(str(directory / "independent_image_planes.png"), dpi=180)
    if result["imaging_parameters"].detector.enabled and detector_figure is not None:
        detector_figure.savefig(str(directory / "detector_pixels.png"), dpi=180)
    arrays = {}
    summaries = []
    for index, channel in enumerate(result["channels"], 1):
        summary = {key: value for key, value in channel.items()
                   if not isinstance(value, (np.ndarray, dict))}
        if "detector" in channel:
            summary["detector"] = {key: value for key, value in channel["detector"].items()
                                   if not isinstance(value, np.ndarray)}
            arrays.update({"ch{}_detector_{}".format(index, key): value
                           for key, value in channel["detector"].items()
                           if isinstance(value, np.ndarray)})
        summaries.append(summary)
        arrays.update({"ch{}_{}".format(index, key): value for key, value in channel.items()
                       if isinstance(value, np.ndarray)})
    np.savez_compressed(str(directory / "coherent_imaging.npz"), **arrays)
    if result['imaging_parameters'].detector.enabled:
        from camera_readout import export_raw_camera
        export_raw_camera(directory, [channel['detector'] for channel in result['channels']])
    payload = {key: asdict(result[key]) for key in
               ("parameters", "objective_parameters", "imaging_parameters")}
    payload.update(channels=summaries, tube_focal_length_mm=result["tube_focal_length_mm"],
                   model_note=result["model_note"], thickness_note=result["thickness_note"],
                   coordinate_convention=incidence_metadata(result["parameters"].incidence_deg),
                   assumptions="Ideal aplanatic air objectives; single sphere; uniform "
                   "plane-wave illumination; particle-only signal; no background leakage, "
                   "roughness, or full near-field solution. Optional detector uses area-integrated "
                   "sampled irradiance, Poisson photo/dark/background counts, full-well clipping "
                   "before Gaussian read noise, then linear gain, offset, ADC rounding/clipping. "
                   "Optional central-square active area and fixed QE/dark maps. No blooming or pixel crosstalk. "
                   "Detector SNR is pre-saturation with known mean background, not a measured SNR. "
                   "E and H*Z0 modes are flux-normalized, not SI fields; s,p inputs "
                   "are incoherent for unpolarized light. FFT images are periodic.")
    (directory / "imaging_parameters.json").write_text(json.dumps(payload, indent=2), encoding="utf-8")
