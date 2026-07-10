"""Compatibility exports for far-field and objective-collection backends."""

from .native_far_field import SCATMECH_IMPORT_ERROR, calculate_far_field
from .objective_collection import calculate_objective_collection_example

__all__ = [
    "SCATMECH_IMPORT_ERROR",
    "calculate_far_field",
    "calculate_objective_collection_example",
]
