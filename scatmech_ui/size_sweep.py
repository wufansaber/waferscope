"""Particle-size sweep calculation and export helpers."""

import csv
import json
from dataclasses import asdict, replace

import numpy as np

from .coordinate_geometry import incidence_metadata
from .native_far_field import calculate_far_field
from .objective_math import objective_collection_masks, sample_far_field


def calculate_particle_size_sweep(
    simulation_parameters,
    objective_parameters,
    sweep_parameters,
    progress_callback=None,
):
    simulation_parameters.validate()
    objective_parameters.validate()
    sweep_parameters.validate()

    sampling = _collection_sampling_grid(objective_parameters)
    diameters = np.asarray(sweep_parameters.diameters_nm, dtype=float)
    channel_values = np.zeros((len(diameters), 3), dtype=float)
    hemisphere_values = np.zeros(len(diameters), dtype=float)
    thickness_note = ""
    model_note = ""

    for index, diameter_nm in enumerate(diameters):
        parameters = replace(
            simulation_parameters,
            particle_diameter_nm=float(diameter_nm),
        )
        far_field = calculate_far_field(
            parameters,
            _diameter_progress_callback(
                progress_callback,
                index,
                len(diameters),
            ),
        )
        sampled_intensity = sample_far_field(
            far_field,
            sampling["kx"],
            sampling["ky"],
            sampling["kz"],
        )["intensity"]
        channel_values[index, :] = [
            _integrate_dsc(sampled_intensity, mask, sampling)
            for mask in sampling["channel_masks"]
        ]
        hemisphere_values[index] = _integrate_dsc(
            sampled_intensity,
            None,
            sampling,
        )
        thickness_note = far_field["thickness_note"]
        model_note = far_field["model_note"]

    if progress_callback is not None:
        progress_callback(100)

    energy_fractions = np.divide(
        channel_values,
        hemisphere_values[:, None],
        out=np.zeros_like(channel_values),
        where=hemisphere_values[:, None] > 0.0,
    )
    return {
        "plot_kind": "particle_size_sweep",
        "parameters": simulation_parameters,
        "simulation_parameters": simulation_parameters,
        "objective_parameters": objective_parameters,
        "sweep_parameters": sweep_parameters,
        "diameters_nm": diameters,
        "channel_collected_dsc_um2": channel_values,
        "hemisphere_dsc_um2": hemisphere_values,
        "channel_energy_fractions": energy_fractions,
        "thickness_note": thickness_note,
        "model_note": model_note,
    }


def _collection_sampling_grid(objective_parameters):
    theta_axis = np.deg2rad(np.linspace(0.0, 89.0, 179))
    phi_axis = np.linspace(0.0, 2.0 * np.pi, 721)
    phi_grid, theta_grid = np.meshgrid(phi_axis, theta_axis)
    sin_theta = np.sin(theta_grid)
    kx = sin_theta * np.cos(phi_grid)
    ky = sin_theta * np.sin(phi_grid)
    kz = np.cos(theta_grid)
    masks = objective_collection_masks(objective_parameters, kx, ky, kz)
    return {
        "theta_axis": theta_axis,
        "phi_axis": phi_axis,
        "sin_theta": sin_theta,
        "kx": kx,
        "ky": ky,
        "kz": kz,
        "channel_masks": masks["channel_masks"],
    }


def _diameter_progress_callback(callback, index, point_count):
    if callback is None:
        return None

    def report(local_progress):
        overall = 100.0 * (index + local_progress / 100.0) / point_count
        callback(int(round(overall)))

    return report


def _integrate_dsc(intensity, mask, sampling):
    values = intensity * sampling["sin_theta"]
    if mask is not None:
        values = np.where(mask, values, 0.0)
    phi_integral = np.trapz(values, sampling["phi_axis"], axis=1)
    return float(np.trapz(phi_integral, sampling["theta_axis"]))


def export_particle_size_sweep(result, directory, figure):
    directory.mkdir(parents=True, exist_ok=True)
    figure.savefig(
        str(directory / "particle_size_sweep.png"),
        dpi=180,
        facecolor="white",
    )
    np.savez_compressed(
        str(directory / "particle_size_sweep_data.npz"),
        diameters_nm=result["diameters_nm"],
        channel_collected_dsc_um2=result["channel_collected_dsc_um2"],
        hemisphere_dsc_um2=result["hemisphere_dsc_um2"],
        channel_energy_fractions=result["channel_energy_fractions"],
    )
    _export_sweep_csv(result, directory / "particle_size_sweep.csv")
    payload = {
        "simulation": asdict(result["simulation_parameters"]),
        "objective": asdict(result["objective_parameters"]),
        "sweep": asdict(result["sweep_parameters"]),
        "coordinate_convention": incidence_metadata(
            result["simulation_parameters"].incidence_deg
        ),
    }
    payload["simulation"]["illuminated_area_um2"] = result[
        "simulation_parameters"
    ].illuminated_area_um2
    payload["wafer_stack_note"] = result.get("thickness_note", "")
    payload["scatmech_model_note"] = result.get("model_note", "")
    with (directory / "particle_size_sweep_parameters.json").open(
        "w",
        encoding="utf-8",
    ) as stream:
        json.dump(payload, stream, indent=2)


def _export_sweep_csv(result, output_path):
    with output_path.open("w", newline="", encoding="utf-8-sig") as stream:
        writer = csv.writer(stream)
        writer.writerow(
            [
                "diameter_nm",
                "ch1_collected_dsc_um2",
                "ch2_collected_dsc_um2",
                "ch3_collected_dsc_um2",
                "upper_hemisphere_dsc_um2",
                "ch1_energy_fraction",
                "ch2_energy_fraction",
                "ch3_energy_fraction",
            ]
        )
        for index, diameter_nm in enumerate(result["diameters_nm"]):
            writer.writerow(
                [
                    diameter_nm,
                    *result["channel_collected_dsc_um2"][index, :],
                    result["hemisphere_dsc_um2"][index],
                    *result["channel_energy_fractions"][index, :],
                ]
            )
