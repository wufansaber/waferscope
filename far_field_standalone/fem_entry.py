"""Frozen solver retaining the source Python -m subprocess protocol."""
import multiprocessing
import os
from pathlib import Path
import sys


def main():
    multiprocessing.freeze_support()
    if not getattr(sys, 'frozen', False):
        sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
    if len(sys.argv) < 3 or sys.argv[1] != '-m':
        raise SystemExit('Use -m fullwave_prototype.guarded_run or a supported FEM module')
    module = sys.argv[2]
    sys.argv = [sys.argv[0]]+sys.argv[3:]
    os.environ.setdefault('MPLBACKEND', 'Agg')
    if module == 'fullwave_prototype.guarded_run':
        from fullwave_prototype.guarded_run import main as run
    elif module == 'fullwave_prototype.imaging_run':
        from fullwave_prototype.imaging_run import main as run
    elif module == 'fullwave_prototype.detection_replay':
        from fullwave_prototype.detection_replay import main as run
    elif module == 'fullwave_prototype.component_compare':
        from fullwave_prototype.component_compare import main as run
    elif module == 'fullwave_prototype.haze_ensemble':
        from fullwave_prototype.haze_ensemble import main as run
    elif module == 'fullwave_prototype.whole_haze':
        from fullwave_prototype.whole_haze import main as run
    elif module == 'fullwave_prototype.hybrid_run':
        from fullwave_prototype.hybrid_run import main as run
    elif module == 'fullwave_prototype.local_spm_validation':
        from fullwave_prototype.local_spm_validation import main as run
    elif module == 'fullwave_prototype.hybrid_range_validation':
        from fullwave_prototype.hybrid_range_validation import main as run
    elif module == 'fullwave_prototype.surface_mesh_audit':
        from fullwave_prototype.surface_mesh_audit import main as run
    elif module == 'fullwave_prototype.illumination_definition_audit':
        from fullwave_prototype.illumination_definition_audit import main as run
    elif module == 'fullwave_prototype.range_refinement':
        from fullwave_prototype.range_refinement import main as run
    else:
        raise SystemExit('Unsupported packaged module: '+module)
    run()


if __name__ == '__main__':
    main()
