# Initial native-Windows validation, 2026-09-25

Subsequent physical audit: [PHYSICAL_AUDIT_CN.md](PHYSICAL_AUDIT_CN.md).
It found a 4.78% free-sphere energy-balance discrepancy at default order 2
and 1.45%-1.61% rough-patch PML-strength sensitivity. Earlier mesh and FFT
agreement numbers below are not overall accuracy bounds.

Environment: NGSolve/Netgen 6.2.2607, Windows x64, Python 3.11.16, two NGSolve
threads. Independent reference: pySCATMECH 0.1.9 in the unchanged NIST2 environment.
All seven cases below were actually run. They are not theoretical estimates.

Physics: 355 nm illumination, 100 nm PSL sphere (n=1.59, k=0), substrate
n=5.57, k=0.387, incidence 75 deg from normal. These are existing application
input values, not a new material calibration. Reference multipole order is 14.

| Case | Gap nm | Pol | FEM order | DOFs | Upper DSC um^2 | Relative discrepancy |
| --- | ---: | --- | ---: | ---: | ---: | ---: |
| Initial mesh | 2 | p | 1 | 24500 | 0.0012151279 | 11.164% |
| Same mesh | 2 | p | 2 | 95415 | 0.0013436095 | 1.771% |
| Contact | 0 | p | 2 | 95748 | 0.0014189466 | 1.774% |
| Contact, mesh sizes x0.8 | 0 | p | 2 | 159078 | 0.0014227606 | 1.510% |
| Contact, PML thickness 0.24 um | 0 | p | 2 | 137469 | 0.0014034134 | 2.849% |
| Contact, s polarization | 0 | s | 2 | 95748 | 0.00019922953 | 2.384% |
| Contact, third order | 0 | p | 3 | 243744 | 0.0014375658 | 0.485% |

The contact p reference integral is 0.00144457182777 um^2. All integrals use
the same 5 degree angular grid, so this table does not quantify angular
quadrature error. The weighted angular DSC discrepancy in the third-order
case is 0.480%; increasing its sphere radiation quadrature from radial order
10 to 14 changes the sampled DSC by a relative L2 of 2.46e-5.

The larger-PML case also regenerated its mesh, so its change cannot be
attributed solely to PML truncation. This is evidence that the total numerical
error budget remains open, not a reason to choose only the best-looking row.

Third-order solve time was about 101 seconds (103 seconds including field
sampling and reference comparison, excluding figure rendering); an observed
Python working set during the solve was about 16 GiB. Default order remains
2 to reduce resource use. No EXE packaging was attempted.

## Checks passed

Eight tests in the native FEM environment cover parameter validation, Fresnel
limits, tangential E/H continuity, substrate absorption direction, sphere
quadrature moments, free-space dipole radiation, translation phase, and a
native zero-contrast solve. The existing NIST2 suite reported 39 tests, 38
passing and one native-FEM test skipped; that skipped test passed separately
in the FEM environment. All 31 pre-existing tests still pass.

## Remaining validation gates

- Separate mesh/order, physical domain and PML convergence; repeat for s polarization.
- Near-field convergence away from interfaces and at the contact region; no
  claim of sub-percent accuracy for local peak field enhancement.
- Full energy/absorption balance and complex far-field phase-reference alignment.
- Material/wavelength, particle-size and incidence-angle coverage.
- Generalized geometry/radiation handling before adding rough surfaces or layers.

## Rough-patch extension, 2026-09-25

The geometry/radiation extension is now implemented for a finite, triangulated
rough patch (not films). See [ROUGHNESS.md](ROUGHNESS.md). The following are
actual paired solves, using the same mesh for each particle-present/absent pair.
Conditions are the optical parameters above, sphere bottom z=10 nm, a 200 x
200 nm patch of RMS 2 nm, 5 x 5 height nodes, smoothing 0.6 cells, seed 7.
The surface ranges from -4.28747 to +6.55153 nm. Both contrast volumes are
2.96573354e-5 um^3, consistent with zero projected-area mean height.

All four columns below are upper-hemisphere integrated DSC in um^2, not power
percentages. Delta means the difference of the paired complex fields.

| Case | DOFs | Angle step | With particle | Without particle | Delta-field power | Interference |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Order 1 | 29738 | 15 deg | 0.00118395069 | 9.02034689e-6 | 0.00109104647 | 8.38838693e-5 |
| Order 2 | 116706 | 5 deg | 0.00137735684 | 9.79815657e-6 | 0.00127554039 | 9.20182884e-5 |
| Order 2, mesh sizes x0.8 | 179427 | 5 deg | 0.00137887487 | 9.72459623e-6 | 0.00127600907 | 9.31412045e-5 |

The order-1 row also has coarser angular sampling, so its difference is not
solely an element-order effect. Between the order-2 meshes, relative to the
refined result, changes are 0.1101% total, 0.7564% background, 0.0367%
delta-field power and 1.2056% interference. This single refinement is NOT an
overall accuracy bound or proof of asymptotic convergence.

On the order-2 base mesh, subsampling the same 5-degree data at 10 and 15
degrees gives total integrals 0.00136712454 and 0.00133908411 um^2.
The 10-to-5-degree change is about 0.74%, exceeding the observed total-power
mesh change. Finer angular integration is still needed.

The refined pair took 147.6 seconds including mesh/source extraction (excluding
plot rendering). A sampled working set during the run was about 9 GiB; this is
not a guaranteed peak. The base order-2 pair took 59.7 seconds. Native volume
quadrature order 5 versus 7 changed sampled complex amplitudes by relative
L2 values 5.58e-8 (without) and 1.87e-8 (with) in the refined pair. This only
checks extraction quadrature, not FEM error. The signed interference is 6.75%
of the total integral for this realization; it is not a collection efficiency.

The smooth-surface runner was rechecked after the extension at gap 10 nm:
FEM order 2 gives 0.00114174474 um^2 versus SCATMECH 0.00116124125 um^2,
a 1.679% discrepancy on the same 5-degree angular grid. This is a smooth
reference check, NOT an independent validation of rough-surface scattering.

Seventeen analytic/native full-wave tests pass, including geometry volume,
triangle RMS, same-material paired cancellation, source-quadrature comparison,
reciprocal interface conditions, the critical-angle limit and signed interference.
The complete NIST2 suite has 48 tests: 44 pass, four native-FEM tests skip;
those four run and pass separately in the native environment.

Outputs are under `build/fem_roughness/order1`, `order2`, `refined` and
`flat_reference`. Roughness-specific remaining gates include another independent
solver/analytic perturbative benchmark, PML/domain convergence, both polarizations,
local near-field convergence, global energy balance, statistically meaningful
large-patch/ensemble behavior, and the three-channel image-plane connection.

## Near/far/image pipeline, 2026-09-26

The three-channel ideal imaging connection is now implemented, without
replacing the production UI. Full-wave near/far fields were recomputed using
the order-2 rough-patch case above, and the saved complex volume sources were
sampled directly on each pupil. No angular-intensity interpolation is used.
NA is 0.5 for all channels; tilts are +45, -45, 0 degrees about global Y.
Magnification is 20, objective focal length 20 mm, tube focal length 400 mm,
incident irradiance 1e-6 W/um^2, transmission 1, no analyzer or extra defocus.

| Channel | 64 x 64 pupil, power pW | 128 x 128 pupil, power pW | Relative change vs 128 |
| --- | ---: | ---: | ---: |
| +45 deg | 172.692327 | 172.394620 | 0.1727% |
| -45 deg | 210.676546 | 210.328206 | 0.1656% |
| 0 deg | 58.090857 | 57.943288 | 0.2547% |

Powers integrate the full FFT image, not just the displayed central +/-15 um
window. Channels are independent and must not be summed as a union collection
efficiency. The 64-sample render uses padding 8 and the 128-sample check uses
padding 2; padding changes image sampling, not integrated power. Thus this is
a pupil-integral convergence check, not a pointwise comparison of image peaks.
FFT edge absolute-flux fractions are 0.42%, 0.40%, 0.62% for 64 samples and
0.20%, 0.20%, 0.30% for 128 samples. These are diagnostics, not bounds on aliasing.

Pupil/image discrete power disagreement is below 4e-14 relative, and coherent
image-decomposition residuals are below 3e-16 relative to the image maximum.
Neither number estimates physical model accuracy. Seven adapter tests cover
NA=0, very small positive NA, analytic dipole-cap collection, inversion under
source translation, irradiance/transmission scaling, analyzers, magnification,
defocus and coherent image decomposition. Existing vector propagation tests
also still pass after moving their implementation into optical_propagation.py.
Full NIST2 regression: 55 tests, 51 pass and four native-FEM tests skip.

Results: build/fullwave_imaging/default and build/fullwave_imaging/pupil128.
Remaining gates: FEM and optical physical benchmarks, PML/domain/near-field
convergence, absorption balance, finer aperture/field sampling, real lens
aberrations, specular leakage, finite illumination, films and detector response.
The end-to-end chain is functional but remains an experimental ideal model.

Generated arrays, PNGs and per-run JSON are under `build/fem_validation/`.
Run `python -m fullwave_prototype.summarize` in the FEM environment to regenerate
the combined `SUMMARY.md` and `summary.json` there. Generated files and the
isolated runtime are ignored by Git; source and this validation record are not.
