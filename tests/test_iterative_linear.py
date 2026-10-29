import numpy as np
import pytest

from fullwave_prototype.iterative_linear import bounded_gmres


def test_complex_system_agrees_with_direct_solution_and_true_residual():
    matrix = np.array([[2+1j, .2], [.2, 3-.4j]])
    rhs = np.array([1+2j, 3-.2j])
    solution, report = bounded_gmres(matrix, rhs, maxsteps=10)
    np.testing.assert_allclose(solution, np.linalg.solve(matrix, rhs), rtol=1e-10)
    assert report['scipy_info'] == 0
    assert report['scipy_true_relative_residual'] < 1e-10
    assert report['iterations'] <= report['effective_maxsteps'] <= 10


def test_restart_cannot_reset_the_total_iteration_budget():
    matrix = np.diag(np.linspace(1., 100., 100))
    _, report = bounded_gmres(matrix, np.ones(100), maxsteps=83)
    assert report['iterations'] <= report['effective_maxsteps'] <= 83
    assert report['effective_maxsteps'] == 80


def test_unconverged_solution_is_not_reported_as_success():
    _, report = bounded_gmres(np.diag(np.arange(1., 101.)), np.ones(100), maxsteps=3)
    assert report['scipy_info'] > 0
    assert report['scipy_true_relative_residual'] > 1e-8


@pytest.mark.parametrize('options', [{'maxsteps': 0}, {'tolerance': float('nan')}, {'tolerance': -1.}])
def test_invalid_solver_limits_are_rejected(options):
    with pytest.raises(ValueError):
        bounded_gmres(np.eye(2), np.ones(2), **options)
