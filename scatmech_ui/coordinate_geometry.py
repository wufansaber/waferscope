"""Coordinate conventions shared by far-field calculations and plots."""

from math import cos, radians, sin


INCIDENCE_PLANE = "KY-KZ"
ILLUMINATION_ROTATION_DEG = 0.0
SCATMECH_AZIMUTH_OFFSET_DEG = 90.0
INCIDENCE_METADATA_KEYS = (
    "illumination_rotation_deg",
    "scatmech_azimuth_offset_deg",
    "incidence_plane",
    "incident_direction",
    "incident_source_direction",
    "incident_projection_direction",
)


def incident_propagation_direction(incidence_deg):
    """Return the unit propagation vector from the -KY source to the wafer."""
    incidence = radians(incidence_deg)
    return (0.0, sin(incidence), -cos(incidence))


def incident_source_direction(incidence_deg):
    """Return the unit vector from the wafer origin toward the light source."""
    propagation = incident_propagation_direction(incidence_deg)
    return tuple(-component for component in propagation)


def incident_projection_direction(incidence_deg):
    """Return the in-plane component of the incident propagation vector."""
    propagation = incident_propagation_direction(incidence_deg)
    return (propagation[0], propagation[1], 0.0)


def incidence_projection_text(incidence_deg):
    """Return a concise label for the -KY illumination convention."""
    if abs(sin(radians(incidence_deg))) < 1e-12:
        return "none (normal incidence)"
    return "+KY; source -KY"


def scatmech_to_global_xy(local_x, local_y):
    """Rotate SCATMECH's incidence-plane axes by +90 degrees about +KZ."""
    return -local_y, local_x


def global_to_scatmech_xy(global_x, global_y):
    """Rotate global XY components back to SCATMECH's local axes."""
    return global_y, -global_x


def incidence_metadata(incidence_deg):
    """Return a serializable record of the -KY illumination convention."""
    return {
        "illumination_rotation_deg": ILLUMINATION_ROTATION_DEG,
        "scatmech_azimuth_offset_deg": SCATMECH_AZIMUTH_OFFSET_DEG,
        "incidence_plane": INCIDENCE_PLANE,
        "incident_direction": incident_propagation_direction(incidence_deg),
        "incident_source_direction": incident_source_direction(incidence_deg),
        "incident_projection_direction": incident_projection_direction(
            incidence_deg
        ),
    }


def incidence_metadata_from_result(result):
    """Extract JSON-safe coordinate metadata from a simulation result."""
    metadata = {}
    for key in INCIDENCE_METADATA_KEYS:
        if key not in result:
            continue
        value = result[key]
        if hasattr(value, "tolist"):
            value = value.tolist()
        metadata[key] = value
    return metadata
