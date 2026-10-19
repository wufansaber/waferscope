"""Standalone entry point for the far-field simulation UI."""

import sys
from pathlib import Path


def _add_project_root_to_path():
    if getattr(sys, "frozen", False):
        return
    project_root = Path(__file__).resolve().parent.parent
    project_root_text = str(project_root)
    if project_root_text not in sys.path:
        sys.path.insert(0, project_root_text)


_add_project_root_to_path()

from scatmech_ui.app import main  # noqa: E402


if __name__ == "__main__":
    main()
