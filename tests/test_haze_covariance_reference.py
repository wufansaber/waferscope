"""Independent dense covariance checks for the finite coherent-beam reference."""
import numpy as np

from fullwave_prototype.parameters import Parameters
from fullwave_prototype.surface_scattering import radiation_kernel
from tools.haze_covariance_reference import covariance_grid, expected_tensor, expected_intensity


def test_fourier_covariance_matches_dense_physical_covariance():
    n, step, sigma, length = 8, .04, .0002, .05
    rng = np.random.default_rng(27)
    fields = rng.normal(size=(n, n, 3))+1j*rng.normal(size=(n, n, 3))
    tensor = expected_tensor(fields, step, sigma, length)
    p = Parameters(wavelength_nm=266)
    xy = np.array([[1, 0], [-1, 1], [0, -1]])/(2*n*step)
    qxy = xy*(p.wavelength_nm/1000)
    q = np.column_stack((qxy, np.sqrt(1-np.sum(qxy**2, axis=1))))
    yy, xx = np.indices((n, n))
    x, y = xx.ravel(), yy.ravel()
    covariance = covariance_grid(2*n, step, sigma, length)[
        (y[:, None]-y[None, :]) % (2*n), (x[:, None]-x[None, :]) % (2*n)]
    expected = []
    for direction, frequency in zip(q, xy):
        phase = np.exp(-2j*np.pi*step*(frequency[0]*x+frequency[1]*y))
        field = phase[:, None]*fields.reshape(-1, 3)
        moment = field.T@covariance@field.conj()*step**4
        kernel = radiation_kernel(direction, p)
        expected.append(np.einsum('aj,jk,ak->', kernel, moment, kernel.conj()).real)
    np.testing.assert_allclose(expected_intensity(tensor, step, q, p), expected, rtol=2e-12)


def test_rms_changes_expected_intensity_quadratically():
    fields = np.ones((8, 8, 3), complex)
    small = expected_tensor(fields, .04, .0002, .05)
    large = expected_tensor(fields, .04, .002, .05)
    for key in small:
        np.testing.assert_allclose(large[key], small[key]*100, rtol=1e-11, atol=1e-24)
