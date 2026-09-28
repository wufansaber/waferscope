# Independent coherent image planes

The three channels have separate objectives, tube lenses, and image planes.
The default objective focal length is 20 mm and magnification is 20, giving
a 400 mm tube focal length. Both sides are air. Objectives obey the sine
condition and are ideal aplanatic systems with unit pupil magnification.
This is a semi-analytic single-particle image model, not a full near-field solver.

## Inputs and workflow

Use the existing particle, wafer stack, illumination, and three NA/tilt controls.
In **Independent image planes**, set magnification, incident irradiance,
image-side defocus, power transmission, analyzer, and pupil samples per axis.
Multipole order 0 selects the SCATMECH automatic estimate; positive values
set lmax explicitly for convergence checks. Repeated sphere/surface scattering
uses order=-1 and Norm_Inc_Approx=0.
Click **Calculate three image planes**. Each run constructs a fresh SCATMECH
model from current inputs; it does not reuse an older intensity plot.
NA=0 disables a channel. Optical settings in this group apply to all three
independent channels. The analyzer axes are local collimated-pupil axes.

Incident irradiance is power per area perpendicular to the incident beam,
not power per wafer footprint area. Output image irradiance is W/um^2 and
channel power is W. No per-channel peak normalization is used. Each image
has its own absolute color scale with a power-law display stretch (gamma 0.35).
The line profiles remain linear. Displayed images crop the central PSF;
exports preserve the complete calculated field of view.

## Field and energy conventions

SCATMECH uses exp(-i omega t). Python JonesDSC rows and columns are ordered
(s,p), with s = unit(z cross k) and p = k cross s. The native incidence plane
is rotated +90 degrees to the application's KY-KZ incidence plane.
For unit incident polarization e, the DSC is norm(J e)^2 in um^2/sr.

For local object-side direction (u,v,w), w=sqrt(1-u^2-v^2), the measure is
dOmega=du dv/w. Complex Cartesian fields are rotated into the channel frame
and parallel-transported onto the collimated pupil. The resulting amplitude
Q includes 1/sqrt(w), so integral norm(Q)^2 du dv is collected DSC.
No extra cos(incidence) multiplies particle DSC.

The tube lens maps to image directions (-u/M,-v/M,sqrt(1-(u^2+v^2)/M^2)).
Polarizations are rotated onto these directions and weighted by the inverse
square root of the image direction cosine to preserve axial flux. Image
fields are evaluated as a discrete vector angular-spectrum integral:

    A(x,y) = 1/(lambda M) integral P(u,v)
             exp[-i 2 pi (u x + v y)/(lambda M)] du dv

Here P includes the rotations, amplitude transmission, and image defocus
phase exp[i k dz (w_image-1)]. H*Z0 is propagated alongside E; image axial
flux is Re(E_x (H*Z0)_y* - E_y (H*Z0)_x*). These stored modes are normalized
to incident irradiance, not SI fields. To get SI E, multiply E modes by
sqrt(2 Z0 I_inc), with I_inc in W/m^2; H modes use sqrt(2 I_inc/Z0).

For unpolarized illumination the two incident basis modes are averaged
in intensity, never added coherently. Channels are not combined in field.
Reported powers belong to separate ideal receivers. Overlapping ideal
cones cannot both receive the same physical rays without an explicit
splitting/obstruction model; a sum is not a union collection efficiency.

## Scope and accuracy

- Single sphere above a planar multilayer substrate, uniform plane wave,
  ideal aplanatic lenses, no vignetting, aberrations, or background leakage.
- Images are particle-scattered signal only. A warning is shown if a cone
  admits the specular direction. This does not implement a physical stop
  or guarantee a dark background in a real instrument.
- Roughness-enabled coherent imaging is rejected: the statistical Mueller
  BRDF does not uniquely determine the spatial phase correlations.
- No full near field, particle interior, finite illumination beam, or
  finite-wafer coherence refinement is added. Optional detector sampling
  and noise are described below.
- Material n,k remain the user-specified values. Changing wavelength does
  not automatically supply measured material dispersion.
- The focal reference is the native SCATMECH amplitude reference. No scan
  across a tilted wafer or field-dependent focusing is simulated here.
- FFT sampling is periodic. Increasing pupil samples enlarges image FOV;
  increasing padding refines image spacing, not angular accuracy.
- The reported flux error tests discrete Parseval normalization only.
  It is not an estimate of model error. Check pupil sampling (32/64/96/128)
  and SCATMECH multipole-order convergence separately. Near grazing collection
  needs special care because of the pupil Jacobian.

## Independent detectors

Enable **Detectors (three independent ROIs)** to evaluate a centered rectangular
pixel array for each channel. Camera settings are shared, but each channel has
its own image, signal and independent random-number stream. NA=0 removes optical
signal; a powered detector still produces dark/background/read noise.
Defaults (illustrative, not a calibrated camera) are 6.5 um pixels, 32 by 32 ROI,
1 ms exposure, QE=0.6, 2 electron RMS read noise, 0.1 electron/s/pixel dark current,
zero extra background and 30000 electron full well. Set QE at the laser wavelength.

For each pixel j, sampled image irradiance is integrated over its area:

    P_j = integral_pixel I(x,y) dx dy
    photons_j = P_j t / (h c / lambda)
    S_j = QE photons_j
    B_j = t (dark_current + background_rate)
    charge_j ~ Poisson(S_j + B_j)
    readout_j = min(charge_j, full_well) + Normal(0, read_noise^2)

Integration uses exact cell/pixel overlap weights for the piecewise-constant
sampled image, not a point sample at the pixel center. This does not remove
optical sampling error: refine the image grid/pupil and check convergence.
The detector ROI must fit entirely within the computed image field. No FFT
wrapping or zero extrapolation is allowed. Increase pupil samples to enlarge
the field; padding only refines sample spacing. Materially negative axial flux
is rejected rather than silently interpreted as a photon rate.

Reported pixel and ROI SNR are ideal pre-saturation quantities with a known
mean background (no additional reference-frame estimation noise):

    SNR_j = S_j / sqrt(S_j + B_j + read_noise^2)
    SNR_ROI = sum(S_j) / sqrt(sum(S_j+B_j) + N_pixels read_noise^2)

They are not valid descriptions of a saturated measurement. Both expected
overfull pixels and sampled saturated pixels are exported. Readout can be
negative or exceed full well because read noise is added after charge clipping.
No ADC quantization, electronic offset/gain, blooming, pixel crosstalk, or
wavelength-dependent detector calibration is simulated. Background is an
independent uniform electron rate, not a computed roughness or leakage field.
Signal views and readout views have separate absolute electron-per-pixel scales
(power-law display stretch gamma 0.5), without peak normalization.

## Exports

Export saves a PNG, an NPZ containing the complex Jones pupils and image
E/H modes for both incident polarizations, and JSON with the full parameter
snapshot, channel powers, flux errors, and image-edge fractions. Array axes
are [y,x,Cartesian component,incident polarization]. Jones arrays instead
end in [outgoing s/p,incoming s/p]. Image coordinates are in um; pupil axes
are in mm. Axis origin is the central image sample.

With detectors enabled, the NPZ also contains `ch1_detector_*` through
`ch3_detector_*` arrays: integrated pixel powers, expected photons and electrons,
sampled readout, saturation masks, ideal pixel SNR and pixel axes. Detector arrays
are [row/y,column/x]. JSON contains the detector settings and per-channel summaries.
UI export additionally writes `detector_pixels.png`.

## Sources

- https://pages.nist.gov/pySCATMECH/local.html
- https://pages.nist.gov/SCATMECH/docs/bobvlieg.htm
- https://pages.nist.gov/SCATMECH/code/jmatrix.cpp
- https://doi.org/10.1098/rspa.1959.0200
