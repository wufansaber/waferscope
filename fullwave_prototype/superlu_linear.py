"""Validate full/triangular CSR against the native complex operator before LU."""
import numpy as np
from scipy.sparse import diags


def validated_matrix(csr, action, tolerance=1e-12):
    if csr.shape[0] != csr.shape[1] or not csr.shape[0] or not np.isfinite(csr.data).all():
        raise ValueError('A finite nonempty square matrix is required')
    rng = np.random.default_rng(20261027)
    probes = [rng.normal(size=csr.shape[0])+1j*rng.normal(size=csr.shape[0]) for _ in range(2)]
    reference = [action(x) for x in probes]

    def error(matrix):
        values = [np.linalg.norm(matrix@x-y)/max(np.linalg.norm(y), 1e-300)
                  for x, y in zip(probes, reference)]
        return float(max(values))

    change = error(csr)
    if np.isfinite(change) and change < tolerance:
        return csr.tocsc(), dict(storage='full CSR', operator_relative_l2=change)
    # Maxwell/PML is complex symmetric, not Hermitian: never conjugate here.
    full = csr+csr.T-diags(csr.diagonal())
    change = error(full)
    if not np.isfinite(change) or change >= tolerance:
        raise ValueError('CSR reconstruction does not match the native operator')
    return full.tocsc(), dict(storage='complex-symmetric triangular CSR', operator_relative_l2=change)
