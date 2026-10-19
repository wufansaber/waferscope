"""Band-limited experimental height/PSD inputs with explicit SI-scaled metadata."""
from functools import lru_cache
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.interpolate import RegularGridInterpolator


def input_fingerprint(path):
    source = Path(path).resolve()
    digest = hashlib.sha256(source.read_bytes())
    sidecar = source.with_suffix('.json')
    if sidecar.is_file():
        digest.update(sidecar.read_bytes())
    return digest.hexdigest()


def _metadata(source):
    sidecar = source.with_suffix('.json')
    value = json.loads(sidecar.read_text(encoding='utf-8')) if sidecar.is_file() else {}
    if not isinstance(value, dict):
        raise ValueError('Measurement metadata must be a JSON object')
    return value


def _positive(value, name):
    value = float(value)
    if not np.isfinite(value) or value <= 0:
        raise ValueError(name+' must be positive and finite')
    return value


def load_height(path, detrend='plane', missing='reject'):
    """Rows are Y, columns X; heights are nm, pixel spacing is um."""
    source = Path(path).resolve()
    metadata = _metadata(source)
    if source.suffix.lower() == '.npz':
        with np.load(source, allow_pickle=False) as data:
            if not {'height_nm', 'dx_um', 'dy_um'}.issubset(data.files):
                raise ValueError('Height NPZ requires height_nm, dx_um and dy_um')
            height = np.asarray(data['height_nm'], dtype=float)
            dx, dy = _positive(data['dx_um'], 'dx_um'), _positive(data['dy_um'], 'dy_um')
    elif source.suffix.lower() == '.csv':
        height = np.loadtxt(source, delimiter=',')
        if not {'dx_um', 'dy_um'}.issubset(metadata):
            raise ValueError('Height CSV needs dx_um and dy_um in its same-name JSON sidecar')
        dx, dy = _positive(metadata['dx_um'], 'dx_um'), _positive(metadata['dy_um'], 'dy_um')
    else:
        raise ValueError('Height input must be NPZ or CSV with a JSON pixel-size sidecar')
    if height.ndim != 2 or min(height.shape) < 4 or height.size > 1024**2:
        raise ValueError('Height input needs a 2D grid, at least 4 x 4, at most 1024 x 1024')
    if np.isinf(height).any():
        raise ValueError('Infinite height values are invalid')
    valid = np.isfinite(height)
    if missing not in ('reject', 'nearest') or detrend not in ('none', 'mean', 'plane'):
        raise ValueError('Unknown missing-point or detrending policy')
    absent = int((~valid).sum())
    if absent:
        if missing == 'reject':
            raise ValueError('Height input contains missing points; select an explicit filling policy')
        if valid.sum() < 4 or absent/height.size > .05:
            raise ValueError('Nearest filling requires at least four valid points and <= 5% missing')
        from scipy.ndimage import distance_transform_edt
        indices = distance_transform_edt(~valid, return_distances=False, return_indices=True)
        height = height[tuple(indices)]
    ny, nx = height.shape
    x, y = np.meshgrid(np.arange(nx)*dx, np.arange(ny)*dy)
    coefficients = None
    if detrend == 'mean':
        coefficients = [float(height.mean())]
        height = height-height.mean()
    elif detrend == 'plane':
        design = np.column_stack((np.ones(height.size), x.ravel(), y.ravel()))
        coefficients = np.linalg.lstsq(design, height.ravel(), rcond=None)[0]
        height = height-(design @ coefficients).reshape(height.shape)
        coefficients = coefficients.tolist()
    report = dict(source=str(source), sha256=input_fingerprint(source),
        kind='height', shape_yx=[ny, nx], dx_um=dx, dy_um=dy,
        footprint_xy_um=[(nx-1)*dx, (ny-1)*dy],
        frequency_resolution_xy_per_um=[1/(nx*dx), 1/(ny*dy)],
        nyquist_xy_per_um=[1/(2*dx), 1/(2*dy)],
        height_units='nm', lateral_units='um', detrend=detrend,
        detrend_coefficients=coefficients, missing_policy=missing, missing_points=absent,
        rms_nm=float(np.sqrt(np.mean(height**2))), metadata=metadata,
        instrument_transfer_correction='not applied; source calibration must be supplied externally')
    return height, report


def height_periodogram(height_nm, dx_um, dy_um):
    """Two-sided 2D density; its DFT-bin integral obeys Parseval without a window."""
    height = np.asarray(height_nm, dtype=float)/1000
    if height.ndim != 2 or not np.isfinite(height).all():
        raise ValueError('Periodogram height must be a finite 2D array')
    ny, nx = height.shape
    density = abs(np.fft.fftshift(np.fft.fft2(height)))**2*dx_um*dy_um/(nx*ny)
    return (np.fft.fftshift(np.fft.fftfreq(nx, dx_um)),
            np.fft.fftshift(np.fft.fftfreq(ny, dy_um)), density)


def load_psd(path):
    """NPZ stores fx_per_um, fy_per_um, psd_um4 (Y,X); CSV stores radial f,S2."""
    source = Path(path).resolve()
    metadata = _metadata(source)
    if source.suffix.lower() == '.csv':
        table = np.loadtxt(source, delimiter=',', comments='#', ndmin=2)
        if table.shape[1] != 2 or not np.isfinite(table).all():
            raise ValueError('Radial PSD CSV needs two finite columns: f_per_um, S2_um4')
        frequency, radial = table.T
        if (len(frequency) < 4 or frequency[0] < 0 or np.any(np.diff(frequency) <= 0)
                or np.any(radial < 0) or not np.any(radial > 0)):
            raise ValueError('PSD frequencies must be nonnegative and increasing; density must be nonnegative')
        # Midpoint bins end inside the measured radial bandwidth; never extrapolate.
        count = min(129, max(33, 2*len(frequency)-1))
        if count % 2 == 0:
            count -= 1
        spacing = 2*frequency[-1]/count
        fx = fy = (np.arange(count)-count//2)*spacing
        xx, yy = np.meshgrid(fx, fy)
        density = np.interp(np.hypot(xx, yy), frequency, radial, left=0, right=0)
        radial_data = (frequency, radial)
    elif source.suffix.lower() == '.npz':
        with np.load(source, allow_pickle=False) as data:
            if not {'fx_per_um', 'fy_per_um', 'psd_um4'}.issubset(data.files):
                raise ValueError('PSD NPZ requires fx_per_um, fy_per_um and psd_um4')
            fx, fy = np.asarray(data['fx_per_um'], float), np.asarray(data['fy_per_um'], float)
            density = np.asarray(data['psd_um4'], float)
        radial_data = None
    else:
        raise ValueError('PSD input must be NPZ or a two-column radial CSV')
    if (fx.ndim != 1 or fy.ndim != 1 or min(len(fx), len(fy)) < 4 or
            max(len(fx), len(fy)) > 1024 or density.shape != (len(fy), len(fx)) or
            not all(np.isfinite(v).all() for v in (fx, fy, density)) or
            np.any(density < 0) or not np.any(density > 0)):
        raise ValueError('Invalid two-sided PSD grid (4..1024 bins/axis, nonnegative um4 density)')
    steps = []
    zero = []
    for axis in (fx, fy):
        delta = np.diff(axis)
        if np.any(delta <= 0) or not np.allclose(delta, delta[0], rtol=1e-8, atol=1e-12):
            raise ValueError('PSD frequency axes must be increasing and uniform')
        center = int(np.argmin(abs(axis)))
        if center != len(axis)//2 or not np.isclose(axis[center], 0, atol=1e-12):
            raise ValueError('PSD frequency axes must be centered DFT bins including zero')
        steps.append(float(delta[0]))
        zero.append(center)
    opposite_x = (2*zero[0]-np.arange(len(fx))) % len(fx)
    opposite_y = (2*zero[1]-np.arange(len(fy))) % len(fy)
    if not np.allclose(density, density[np.ix_(opposite_y, opposite_x)], rtol=1e-5,
                       atol=density.max()*1e-10):
        raise ValueError('A real height field requires S(-fx,-fy)=S(fx,fy)')
    variance = float(density.sum()*steps[0]*steps[1])
    report = dict(source=str(source), sha256=input_fingerprint(source), kind='psd',
        psd_units='um4', frequency_units='cycles/um', shape_yx=list(density.shape),
        frequency_steps_xy_per_um=steps, rms_nm=1000*np.sqrt(variance),
        frequency_bounds_xy_per_um=[[float(fx[0]), float(fx[-1])],
                                   [float(fy[0]), float(fy[-1])]],
        out_of_band='zero; no extrapolation', density_integration='DFT bin sum',
        radial_csv=radial_data is not None, metadata=metadata,
        instrument_transfer_correction='not applied; source calibration must be supplied externally')
    if radial_data is not None:
        report['measured_radial_band_per_um'] = [float(radial_data[0][0]), float(radial_data[0][-1])]
    return fx, fy, density, report


def sample_psd(path, frequency):
    source = Path(path)
    points = np.asarray(frequency, dtype=float)
    if source.suffix.lower() == '.csv':
        table = np.loadtxt(source, delimiter=',', comments='#', ndmin=2)
        return np.interp(np.linalg.norm(points, axis=-1), table[:, 0], table[:, 1],
                         left=0, right=0)
    fx, fy, density, _ = load_psd(path)
    interpolate = RegularGridInterpolator((fy, fx), density, bounds_error=False, fill_value=0)
    return interpolate(points[..., ::-1])


@lru_cache(maxsize=4)
def _psd_kernel(path, fingerprint, count, width_um):
    fx, fy, density, metadata = load_psd(path)
    offsets = np.arange(1-count, count)*width_um/(count-1)
    ex = np.exp(2j*np.pi*offsets[:, None]*fx[None, :])
    ey = np.exp(2j*np.pi*offsets[:, None]*fy[None, :])
    df_x, df_y = metadata['frequency_steps_xy_per_um']
    return (ex @ density.T @ ey.T).real*df_x*df_y


def psd_covariance_kernel(roughness):
    return _psd_kernel(roughness.measurement_path, input_fingerprint(roughness.measurement_path),
                       roughness.grid_points, roughness.width_um)


@lru_cache(maxsize=4)
def _psd_factor(path, fingerprint, count, width_um):
    kernel = _psd_kernel(path, fingerprint, count, width_um)
    ix, iy = np.meshgrid(np.arange(count), np.arange(count), indexing='ij')
    ix, iy = ix.ravel(), iy.ravel()
    covariance = kernel[ix[:, None]-ix[None, :]+count-1,
                        iy[:, None]-iy[None, :]+count-1]
    return np.linalg.cholesky(covariance+np.eye(count**2)*kernel[count-1, count-1]*1e-12)


def measured_nodes(roughness):
    if roughness.model == 'measured_psd':
        factor = _psd_factor(roughness.measurement_path, input_fingerprint(roughness.measurement_path),
                             roughness.grid_points, roughness.width_um)
        rng = np.random.default_rng(roughness.seed)
        return (factor @ rng.standard_normal(roughness.grid_points**2)).reshape(
            (roughness.grid_points,)*2)
    height, metadata = load_height(roughness.measurement_path, roughness.detrend,
                                   roughness.missing_policy)
    ny, nx = height.shape
    dx, dy = metadata['dx_um'], metadata['dy_um']
    x = (np.arange(nx)-(nx-1)/2)*dx
    y = (np.arange(ny)-(ny-1)/2)*dy
    if roughness.width_um > min(x[-1]-x[0], y[-1]-y[0])+1e-12:
        raise ValueError('FEM square patch exceeds the measured height footprint; reduce patch width')
    axis = np.linspace(-roughness.width_um/2, roughness.width_um/2, roughness.grid_points)
    xx, yy = np.meshgrid(axis, axis, indexing='ij')
    interpolate = RegularGridInterpolator((y, x), height/1000, bounds_error=True)
    # The footprint check above allows floating-point endpoint roundoff, not extrapolation.
    return interpolate(np.stack((np.clip(yy, y[0], y[-1]),
                                 np.clip(xx, x[0], x[-1])), axis=-1))


def measurement_report(roughness):
    if roughness.model == 'measured_psd':
        fx, fy, density, report = load_psd(roughness.measurement_path)
        nx = 1/(2*roughness.width_um/(roughness.grid_points-1))
        mask = (abs(fy[:, None]) <= nx) & (abs(fx[None, :]) <= nx)
        report['height_grid_band_variance_fraction'] = float(density[mask].sum()/density.sum())
        report['surface_assumption'] = 'Gaussian random height with the measured two-sided PSD'
        return report
    _, report = load_height(roughness.measurement_path, roughness.detrend, roughness.missing_policy)
    report['fem_crop'] = 'Centered square, bilinear node interpolation, then explicit edge taper'
    return report
