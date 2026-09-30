"""Locate reusable scattering for a displayed or relocated result."""
import json
from pathlib import Path


def reusable_source(directory):
    directory = Path(directory)
    candidates = [directory/'scattering', directory]
    try:
        report = json.loads((directory/'imaging_report.json').read_text(encoding='utf-8'))
        recorded = report.get('source_directory') if isinstance(report, dict) else None
        if recorded:
            recorded = Path(recorded)
            candidates.append(recorded if recorded.is_absolute() else directory/recorded)
    except (OSError, ValueError, TypeError):
        pass
    for candidate in candidates:
        try:
            source_report = json.loads((candidate/'report.json').read_text(encoding='utf-8'))
            if (isinstance(source_report, dict) and source_report.get('source_format_version') == 1
                    and (candidate/'paired_fields.npz').is_file()):
                return candidate.resolve()
        except (OSError, ValueError, TypeError):
            continue
    return None
