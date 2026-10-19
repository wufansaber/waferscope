"""SCATMECH wafer-stack, particle, and microroughness model construction."""

from dataclasses import dataclass

import numpy as np

try:
    from pySCATMECH.brdf import BRDF_Model
    from pySCATMECH.fresnel import Film, FilmStack, OpticalFunction
    from pySCATMECH.local import Local_BRDF_Model
    from pySCATMECH.mueller import Polarization

    SCATMECH_IMPORT_ERROR = None
except ModuleNotFoundError as error:
    if error.name not in ("pySCATMECH", "SCATPY"):
        raise
    BRDF_Model = Film = FilmStack = OpticalFunction = None
    Local_BRDF_Model = Polarization = None
    SCATMECH_IMPORT_ERROR = error


@dataclass(frozen=True)
class ScatteringModelBundle:
    model: object
    output_kind: str
    effective_area_um2: float
    thickness_note: str
    model_note: str
    stack_layer_count: int


def abc_psd_a_um4(rms_nm, correlation_length_um, exponent):
    """Return ABC PSD parameter A for the requested integrated RMS height."""
    rms_um = rms_nm / 1000.0
    return (
        rms_um**2
        * correlation_length_um**2
        * (exponent - 2.0)
        / (2.0 * np.pi)
    )


def build_scattering_model(parameters):
    """Build a particle-only DSC model or a particle-plus-roughness BRDF model."""
    if SCATMECH_IMPORT_ERROR is not None:
        raise RuntimeError(
            "pySCATMECH is unavailable in the active Python interpreter."
        ) from SCATMECH_IMPORT_ERROR

    parameters.validate()
    wavelength_um = parameters.wavelength_nm / 1000.0
    substrate, stack, thickness_note, stack_layer_count = _build_wafer_stack(
        parameters,
        wavelength_um,
    )
    particle_parameters = {
        "lambda": wavelength_um,
        "radius": parameters.particle_diameter_nm / 2000.0,
        "sphere": complex(parameters.particle_n, parameters.particle_k),
        "substrate": substrate,
        "type": 0,
    }
    if stack is not None:
        particle_parameters["stack"] = stack

    if not parameters.roughness_enabled:
        model = Local_BRDF_Model(
            "Bobbert_Vlieger_BRDF_Model",
            particle_parameters,
        )
        return ScatteringModelBundle(
            model=model,
            output_kind="dsc",
            effective_area_um2=1.0,
            thickness_note=thickness_note,
            model_note=(
                "SCATMECH Bobbert_Vlieger_BRDF_Model; surface roughness disabled."
            ),
            stack_layer_count=stack_layer_count,
        )

    effective_area_um2 = parameters.illuminated_area_um2
    particle_parameters["density"] = 1.0 / effective_area_um2
    roughness_model, roughness_name = _build_roughness_model(
        parameters,
        wavelength_um,
        substrate,
        stack,
        stack_layer_count,
    )
    validity_note = _roughness_validity_note(parameters, wavelength_um)
    psd_note = ("Gaussian PSD RMS={:.6g} nm, L={:.6g} um".format(
        parameters.roughness_rms_nm, parameters.roughness_correlation_length_um)
        if parameters.roughness_psd == "gaussian" else
        "ABC PSD RMS={:.6g} nm, B={:.6g} um, C={:.6g}".format(
            parameters.roughness_rms_nm, parameters.roughness_correlation_length_um,
            parameters.roughness_exponent))
    if parameters.haze_only:
        return ScatteringModelBundle(
            model=roughness_model, output_kind="brdf_to_dsc",
            effective_area_um2=effective_area_um2, thickness_note=thickness_note,
            model_note=("SCATMECH {} only (no particle); {}; rough area={:.6g} x {:.6g} um "
                        "({:.6g} um^2).{}".format(
                            roughness_name, psd_note, parameters.illumination_length_um,
                            parameters.illumination_width_um, effective_area_um2, validity_note)),
            stack_layer_count=stack_layer_count)
    particle_model = Local_BRDF_Model(
        "Bobbert_Vlieger_BRDF_Model",
        particle_parameters,
    )
    combined_model = BRDF_Model(
        "Two_Source_BRDF_Model",
        {
            "lambda": wavelength_um,
            "substrate": substrate,
            "type": 0,
            "factor1": 1.0,
            "source1": particle_model,
            "factor2": 1.0,
            "source2": roughness_model,
            "correlation": 0.0,
        },
    )
    model_note = (
        "SCATMECH Two_Source_BRDF_Model: Bobbert_Vlieger_BRDF_Model + {} "
        "(incoherent Mueller BRDF sum); {}; rectangular footprint={:.6g} x {:.6g} um "
        "(area {:.6g} um^2).{}"
    ).format(
        roughness_name,
        psd_note,
        parameters.illumination_length_um,
        parameters.illumination_width_um,
        effective_area_um2,
        validity_note,
    )
    return ScatteringModelBundle(
        model=combined_model,
        output_kind="brdf_to_dsc",
        effective_area_um2=effective_area_um2,
        thickness_note=thickness_note,
        model_note=model_note,
        stack_layer_count=stack_layer_count,
    )


def _build_wafer_stack(parameters, wavelength_um):
    wafer_index = complex(parameters.wafer_n, parameters.wafer_k)
    thickness_um = parameters.wafer_thickness_um
    optical_depth = (
        4.0 * np.pi * parameters.wafer_k * thickness_um / wavelength_um
    )
    coherent_limit_um = 100.0 * wavelength_um
    finite_wafer = thickness_um <= coherent_limit_um and optical_depth < 60.0

    films = []
    if finite_wafer:
        substrate = 1.0 + 0.0j
        films.append(
            Film(OpticalFunction(wafer_index), thickness=thickness_um)
        )
        base_note = (
            "Finite coherent wafer layer: {:.6g} um on air backing."
        ).format(thickness_um)
    else:
        substrate = wafer_index
        base_note = (
            "Wafer is optically/geometrically thick at this wavelength and "
            "is modeled as a semi-infinite substrate."
        )

    for layer in parameters.coating_layers:
        films.append(
            Film(
                OpticalFunction(complex(layer.n, layer.k)),
                thickness=layer.thickness_nm / 1000.0,
            )
        )

    if parameters.coating_layers:
        layer_description = ", ".join(
            "{} {:.6g} nm".format(layer.material, layer.thickness_nm)
            for layer in parameters.coating_layers
        )
        coating_note = (
            " {} coating layer(s), ordered wafer-to-surface: {}."
        ).format(len(parameters.coating_layers), layer_description)
    else:
        coating_note = " No coating layers."

    stack = FilmStack(films) if films else None
    return substrate, stack, base_note + coating_note, len(films)


def _build_roughness_model(
    parameters,
    wavelength_um,
    substrate,
    stack,
    stack_layer_count,
):
    if parameters.roughness_psd == "gaussian":
        psd = {None: "Gaussian_PSD_Function",
               "sigma": parameters.roughness_rms_nm/1000.0,
               "length": parameters.roughness_correlation_length_um}
    else:
        psd = {
            None: "ABC_PSD_Function",
            "A": abc_psd_a_um4(
                parameters.roughness_rms_nm,
                parameters.roughness_correlation_length_um,
                parameters.roughness_exponent,
            ),
            "B": parameters.roughness_correlation_length_um,
            "C": parameters.roughness_exponent,
        }
    model_parameters = {
        "lambda": wavelength_um,
        "substrate": substrate,
        "type": 0,
        "psd": psd,
    }
    if stack is None:
        return (
            BRDF_Model("Microroughness_BRDF_Model", model_parameters),
            "Microroughness_BRDF_Model",
        )

    model_parameters["stack"] = stack
    if parameters.roughness_interface_mode == "top_interface":
        model_parameters["this_layer"] = stack_layer_count
        model_name = "Roughness_Stack_BRDF_Model"
    elif parameters.roughness_interface_mode == "correlated_interfaces":
        model_name = "Correlated_Roughness_Stack_BRDF_Model"
    else:
        model_name = "Uncorrelated_Roughness_Stack_BRDF_Model"
    return BRDF_Model(model_name, model_parameters), model_name


def _roughness_validity_note(parameters, wavelength_um):
    rms_um = parameters.roughness_rms_nm / 1000.0
    smoothness_parameter = (
        4.0
        * np.pi
        * rms_um
        * abs(np.cos(np.deg2rad(parameters.incidence_deg)))
        / wavelength_um
    )
    if smoothness_parameter <= 0.3:
        return ""
    return (
        " Warning: the Rayleigh-Rice smooth-surface parameter is {:.3g}; "
        "verify the perturbation-model validity."
    ).format(smoothness_parameter)
