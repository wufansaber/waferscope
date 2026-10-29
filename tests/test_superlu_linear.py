"""A complex-symmetric PML matrix must never be Hermitian-symmetrized."""
import numpy as np
import pytest
from scipy.sparse import csr_matrix, tril, triu
from scipy.sparse.linalg import splu

from fullwave_prototype.superlu_linear import validated_matrix


@pytest.mark.parametrize('storage', [csr_matrix, tril, triu])
def test_native_operator_matches_full_or_triangular_complex_symmetric_storage(storage):
    actual = csr_matrix(np.array([[2+1j, 1-.3j], [1-.3j, -3+.5j]]))
    matrix, report = validated_matrix(storage(actual).tocsr(), lambda x: actual@x)
    assert report['operator_relative_l2'] < 1e-14
    np.testing.assert_allclose(matrix.toarray(), actual.toarray())
    rhs = np.array([1+2j, 2-1j])
    solution = splu(matrix, permc_spec='MMD_AT_PLUS_A', options={'SymmetricMode': True}).solve(rhs)
    np.testing.assert_allclose(actual@solution, rhs, atol=1e-14)


def test_bad_matrix_is_rejected_instead_of_renormalized():
    actual = csr_matrix(np.array([[2+1j, 1-.3j], [1-.3j, -3+.5j]]))
    with pytest.raises(ValueError, match='native operator'):
        validated_matrix(2*actual, lambda x: actual@x)
