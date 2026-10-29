"""Strictly bounded library GMRES with an explicit true-residual report."""
import inspect
import numpy as np
from scipy.sparse.linalg import gmres


def bounded_gmres(matrix, rhs, preconditioner=None, maxsteps=300, tolerance=1e-10, progress=None):
    if not isinstance(maxsteps, int) or not 1 <= maxsteps <= 1000:
        raise ValueError('Use 1..1000 iteration steps')
    if not np.isfinite(tolerance) or not 0 < tolerance < 1:
        raise ValueError('Tolerance must be finite and in (0,1)')
    rhs = np.asarray(rhs, dtype=complex)
    if rhs.ndim != 1 or not np.isfinite(rhs).all():
        raise ValueError('A finite vector RHS is required')
    restart = min(40, maxsteps)
    cycles = maxsteps//restart
    history = []
    def callback(value):
        history.append(float(value))
        if progress:
            progress(len(history), float(value))
    tolerance_key = 'rtol' if 'rtol' in inspect.signature(gmres).parameters else 'tol'
    solution, info = gmres(matrix, rhs, M=preconditioner, restart=restart, maxiter=cycles,
                           atol=0., callback=callback, callback_type='pr_norm', **{tolerance_key: tolerance})
    residual = np.linalg.norm(matrix@solution-rhs)/max(np.linalg.norm(rhs), 1e-300)
    return solution, dict(iterations=len(history), scipy_info=int(info),
        requested_maxsteps=maxsteps, effective_maxsteps=cycles*restart,
        preconditioned_residuals=history, scipy_true_relative_residual=float(residual))
