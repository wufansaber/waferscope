"""Power-covered stationary rough surfaces, separate from local FEM geometry."""
import numpy as np

from .rough_statistics import abc_amplitude, nyquist_variance_fraction


def window(x, y, width, taper):
    def edge(axis):
        distance = np.maximum(0, width/2-np.abs(axis))
        return np.sin(np.pi/2*np.minimum(1, distance/taper))**2
    return edge(y)[:, None]*edge(x)[None, :]


def plan_surface(p, roughness, coverage=.99, variance=.9, taper_um=.025, max_points=2048):
    if p.beam_profile != 'gaussian' or p.beam_waist_um:
        raise ValueError('Whole-beam haze requires Gaussian footprint illumination')
    if roughness.model not in ('gaussian', 'abc'):
        raise ValueError('Whole-beam haze currently supports Gaussian or ABC PSD')
    if not 0 <= roughness.rms_nm < .05*p.wavelength_nm:
        raise ValueError('First-order whole-beam haze requires RMS / wavelength < 0.05')
    if not .9 <= coverage <= .999 or not .5 <= variance <= .999 or not np.isfinite(taper_um) or taper_um <= 0:
        raise ValueError('Invalid coverage, PSD bandwidth or edge taper target')
    length = roughness.correlation_length_um
    if not np.isfinite(length) or length <= 0 or not np.isfinite(roughness.psd_exponent):
        raise ValueError('Positive finite PSD length required')
    if roughness.model == 'abc' and not 2.01 <= roughness.psd_exponent <= 10:
        raise ValueError('ABC exponent must be in [2.01,10]')
    fmax = (np.sqrt(-np.log(1-variance))/(np.pi*length) if roughness.model == 'gaussian' else
        np.sqrt((1-variance)**(-2/(roughness.psd_exponent-2))-1)/length)
    spacing = min(.04, 1/(2*fmax), p.wavelength_nm/1000/8)
    from .surface_planning import coverage_plan
    actual, plan = coverage_plan(p, spacing, coverage, taper_um, max_points)
    plan.update(target_power_coverage=coverage, target_nyquist_variance_fraction=variance,
        nyquist_variance_fraction=nyquist_variance_fraction(
            roughness.model, plan['spacing_um'], length, roughness.psd_exponent))
    return actual, plan


def random_surface(roughness, plan, seed, return_synthesis_parent=False):
    """Gaussian random height from the absolute PSD, without sample RMS rescaling."""
    n, step = plan['grid_points'], plan['spacing_um']
    parent = 2*n
    fx = np.fft.rfftfreq(parent, step)[None, :]
    fy = np.fft.fftfreq(parent, step)[:, None]
    frequency2 = fx*fx+fy*fy
    length, sigma = roughness.correlation_length_um, roughness.rms_nm/1000
    if roughness.model == 'gaussian':
        spectrum = np.pi*sigma*sigma*length*length*np.exp(-(np.pi*length)**2*frequency2)
    else:
        spectrum = abc_amplitude(roughness.rms_nm, length, roughness.psd_exponent)/(
            1+length*length*frequency2)**(roughness.psd_exponent/2)
    rng = np.random.default_rng(seed)
    transformed = np.fft.rfft2(rng.standard_normal((parent, parent)))
    transformed *= np.sqrt(spectrum)/step
    heights = np.fft.irfft2(transformed, s=(parent, parent))
    if return_synthesis_parent:
        return heights
    # Crop away the synthesis boundary; periodicity is not a physical wafer boundary.
    return heights[n//2:n//2+n, n//2:n//2+n]
