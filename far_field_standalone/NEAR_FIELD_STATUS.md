# Near-field implementation status

The current application calculates the complex scattered far field with
SCATMECH, then propagates independently through ideal collection optics to
image planes and optional detectors. It does NOT yet calculate the full
electromagnetic near field around a sphere touching a layered wafer.

## Native API limitation found during inspection

Inspected source: pySCATMECH 0.1.9 source distribution, `scatmech/bobvlieg2.cpp`,
`Bobbert_Vlieger_BRDF_Model::EField`. The distribution SHA256 is
`0912f6521ed57f1f5caf64e5a5836356ff0a507ae11c11fb347039402d3688d3`.

Its field evaluator returns zero when the sphere-centered radius r is less
than the particle radius a or greater than a+delta. For contact (delta=0),
the supported radial shell degenerates to r=a. Calling this API on a general
XZ/XY grid therefore does not give the surrounding physical near field.
Zero returned outside the supported region must not be plotted as a physical
zero field. The particle interior is also excluded.

## Work still required

1. Expose the solved vector-spherical-wave coefficients without discarding phase.
2. Reconstruct the incident, reflected illumination, particle-scattered and
   substrate-reflected particle fields with consistent coordinates and units.
3. Use the layered-medium angular-spectrum/Green function representation where
   the native local regular-wave expansion is not valid; include evanescent terms.
4. Implement interior fields separately and treat sphere/wafer contact carefully.
5. Validate boundary conditions, loss/energy balances, multipole and integration
   convergence, and the far-distance limit against JonesDSC before enabling UI.

The available Python binary does not by itself supply this complete workflow.
No usable C++ compiler was detected during this implementation pass. Installing
a compiler alone will not fix the mathematical domain restriction.
An inverse FFT of propagating far-field amplitudes cannot restore missing
evanescent information and is not a substitute for this near-field solution.

References:

- [NIST Bobbert-Vlieger model](https://pages.nist.gov/SCATMECH/docs/bobvlieg.htm)
- [NIST source](https://pages.nist.gov/SCATMECH/code/bobvlieg2.cpp)
- [pySCATMECH 0.1.9 source distribution](https://files.pythonhosted.org/packages/fe/ac/43255fc075d48bb0e64bc11c9022515a2254e5242ff50441d06370b72415/pySCATMECH-0.1.9.tar.gz)
