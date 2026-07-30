"""Coordinate conventions shared by far-field calculations and plots."""

from math import cos, radians, sin


INCIDENCE_PLANE = "KX-KZ"
ILLUMINATION_ROTATION_DEG = 90.0
INCIDENCE_METADATA_KEYS = (
    "illumination_rotation_deg",
    "incidence_plane",
    "incident_direction",
    "incident_source_direction",
    "incident_projection_direction",
)


def incident_propagation_direction(incidence_deg):
    """Return the confirmed A-direction unit vector toward the wafer."""
    incidence = radians(incidence_deg)
    return (sin(incidence), 0.0, -cos(incidence))


def incident_source_direction(incidence_deg):
    """Return the unit vector from the wafer origin toward the light source."""
    propagation = incident_propagation_direction(incidence_deg)
    return tuple(-component for component in propagation)


def incident_projection_direction(incidence_deg):
    """Return the in-plane component of the incident propagation vector."""
    propagation = incident_propagation_direction(incidence_deg)
    return (propagation[0], propagation[1], 0.0)


def incidence_projection_text(incidence_deg):
    """Return a concise label for the confirmed illumination convention."""
    if abs(sin(radians(incidence_deg))) < 1e-12:
        return "none (normal incidence)"
    return "+KX; source -KX"


def incidence_metadata(incidence_deg):
    """Return a serializable record of the confirmed A convention."""
    return {
        "illumination_rotation_deg": ILLUMINATION_ROTATION_DEG,
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
