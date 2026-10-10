# Native Windows full-wave prototype

Experimental frequency-domain vector Maxwell solver using Netgen/NGSolve.
No WSL, Linux, commercial solver or change to NIST2 is required. This directory
uses a separate interpreter. The source `far_field_standalone` UI now includes an
experimental FEM tab; the original SCATMECH calculation path is unchanged.
Actual first-run results and remaining accuracy checks: [VALIDATION.md](VALIDATION.md).
Physical audit and failed accuracy gates: [PHYSICAL_AUDIT_CN.md](PHYSICAL_AUDIT_CN.md).
Solver choice and improved numerical controls: [SOLVER_DECISION_CN.md](SOLVER_DECISION_CN.md).
Completed follow-up, source UI and resource limits: [FOLLOWUP_VALIDATION_CN.md](FOLLOWUP_VALIDATION_CN.md).
Equations and assumptions: [THEORETICAL_MODEL_CN.md](THEORETICAL_MODEL_CN.md).
Independent theory checks and balanced settings: [THEORY_VALIDATION_CN.md](THEORY_VALIDATION_CN.md).
Optional lower-order PML and local mesh refinement: [LOW_MEMORY_REFINEMENT_CN.md](LOW_MEMORY_REFINEMENT_CN.md).
The CLI's legacy order-2 default does not establish sub-percent physical accuracy; the
free-sphere energy-balance discrepancy reached 4.78% in the audited case.
The source FEM UI now defaults to the tested `balanced` order-3 configuration,
128 pupil samples and a 12 GiB sampled memory budget. This is not a universal accuracy guarantee.

## Implemented scope

Whole-illumination roughness now has a separate coherent first-order surface
solver and UI entry, with power/taper coverage, PSD sampling, direct aperture
integration and independent image planes. This is NOT full-wafer FEM or
particle-coupled scattering. See [whole-beam haze](WHOLE_BEAM_HAZE_CN.md).

- One homogeneous isotropic sphere in air, touching or separated from a smooth,
  isotropic, absorbing/nonabsorbing semi-infinite substrate at z=0.
- Plane-wave illumination from negative KY toward positive KY and downward Z;
  incidence measured from the surface normal, separate s or p runs.
- Complex total/scattered E and total H*Z0 in a YZ slice, including particle interior.
- Upper-hemisphere scattered amplitudes and differential scattering cross section.
- Independent SCATMECH comparison through the existing NIST2 interpreter.
- A separate finite rough-patch runner with paired particle-present/absent
  solves, coherent interference, and sources above/below the reference plane.
  Gaussian/ABC height statistics and haze ensemble comparisons are described
  in [ROUGHNESS.md](ROUGHNESS.md); research priorities and evidence are in
  [ROUGHNESS_RESEARCH_CN.md](ROUGHNESS_RESEARCH_CN.md).
- A near/far/image pipeline with three independent ideal 20x objectives, direct
  complex FEM pupil sampling and coherent image decomposition. See
  [IMAGING.md](IMAGING.md) for commands, units and limitations.

Finite coherent Gaussian wafer footprints now enter the FEM background, near/far/image
pipeline and haze averaging. The source UI defaults to this mode; the CLI keeps its
plane-wave default. See [coherent illumination](COHERENT_ILLUMINATION_CN.md) for power
normalization, clipping, patch coverage and quantitative limitations. Camera pixel/noise
integration is available in the imaging pipeline, without experimental calibration.

NOT implemented: multilayer FEM films, validated wafer-scale roughness convergence,
other particle shapes or experimentally calibrated camera/illumination response.
Arbitrary-geometry FEM is the longer-term route, not a claim that these
features already work here. Near-field contact peaks are not yet quantitatively
converged. Do not use this prototype as a production replacement.

## Environment

The tested interpreter is a separate Windows Python 3.11 environment, created
under ignored `build/ngsolve_env`. Existing NIST2 remains unchanged.
From the repository root, a fresh setup using an existing Python 3.11 is:

```powershell
& 'D:\anaconda\envs\docling\python.exe' -m venv build\ngsolve_env
& '.\build\ngsolve_env\Scripts\python.exe' -m pip install --only-binary=:all: -r fullwave_prototype\requirements.txt
```

On another machine substitute an installed 64-bit Python 3.11 for the first
path. The current venv depends on that base Python; it is not a portable EXE.
The 2026-09-28 Windows release separately freezes the GUI and solver; it does not
redistribute this nonportable venv. The packaged GUI-to-FEM calculation has been
tested with Anaconda removed from PATH and a working directory outside the repo.
See [packaging verification](../far_field_standalone/PACKAGING_VALIDATION_CN.md).

## Run

```powershell
& '.\build\ngsolve_env\Scripts\python.exe' -m fullwave_prototype.run --output build/fem_validation/default --reference-python D:\anaconda\envs\NIST2\python.exe
```

Or run `fullwave_prototype\run_validation.bat` from a terminal. Without
`--reference-python`, only FEM runs. All outputs remain local. Default physics
matches current application defaults: wavelength 355 nm, sphere diameter 100 nm,
particle n=1.59/k=0, substrate n=5.57/k=0.387, p polarization, incidence 75 deg,
zero gap. These optical constants are parameter values, not a measured dispersion
database. Do not reuse them at 266 nm without checking the materials.

Useful controls: `--gap-nm`, `--polarization s`, `--order 1/2/3`,
`--mesh-scale 0.8` (smaller means finer), `--radius 0.3`,
`--pml-thickness 0.24`, `--wavelength-nm`, `--diameter-nm`, `--incidence-deg`,
`--substrate-n`, `--substrate-k`. Lengths of domain/PML/mesh are in um.
A degree-of-freedom guard stops oversized direct solves; it is not a memory guarantee.

## Mathematical formulation

Time convention is exp(-i omega t); all materials are nonmagnetic.
Let eps_b be air above the interface and substrate permittivity below it.
The analytic background E_b is the incident plane wave plus Fresnel reflection
in air and the transmitted plane wave in the substrate. The contrast delta_eps
is nonzero only within the particle in the smooth benchmark; the rough runner
also includes added substrate above z=0 and removed substrate below z=0.
The FEM scattered field satisfies

    curl curl E_s - k0^2 eps E_s = k0^2 delta_eps E_b
    E_total = E_b + E_s

Use complex H(curl) elements on a curved conforming tetrahedral mesh. A radial
complex-coordinate PML centered on the planar interface preserves its z=0 plane;
outer tangential E_s is set to zero. The PML contains the continued air/substrate
materials, not vacuum on both sides. H_s*Z0 = curl(E_s)/(i k0).

Far fields are reconstructed from the solved polarization delta_eps E_total
inside the sphere, using the air/half-space radiation Green tensor. For an upward
unit direction k+, let k- be its downward counterpart, s the shared azimuthal
basis, and p+/- = k+/- cross s. Define

    F+/- = integral_particle delta_eps E_total(r) exp(-i k0 k+/- dot r) dV
    A_s = k0^2/(4 pi) s dot (F+ + r_s F-)
    A_p = k0^2/(4 pi) (p+ dot F+ + r_p p- dot F-)
    DSC = |A_s|^2 + |A_p|^2

No conjugation is applied in these radiation projections. The reflected term
is essential. This extractor is specific to the flat half-space background and
must be extended before adding films. The rough runner extends the extractor
below z=0 using reciprocal transmitted waves and integrates sources on the
curved FEM elements; its formulas and limitations are in ROUGHNESS.md.
Sphere volume integration uses Gauss radial/polar nodes and uniform azimuthal
nodes, with a second quadrature order check on a subset of directions.

Fields are normalized to unit incident electric amplitude, NOT SI V/m.
DSC is um^2/sr; the hemisphere integral is um^2. At incident irradiance I in
W/um^2, collected power is I times the integrated DSC. Field magnitude squared
in the near field is not equivalent to outward collected power.

## Outputs and validation

The UI also offers a shared-surface coherent hybrid: whole-beam first-order SPM
plus local paired FEM minus the same local triangular SPM amplitudes. See
[hybrid model and validation limits](HYBRID_MODEL_CN.md). It includes same-realization
outer-surface diagnostics and three independent images. Near fields remain local FEM;
it is not a rigorous full-wafer Maxwell solution or an ensemble average.

- `comparison.png`: YZ total/scattered field enhancement and two far-field cuts.
- `fields.npz`: complex field components, outgoing s/p amplitudes, angular DSC
  and coordinate axes. Near arrays are [z,y,component], far arrays [theta,phi,s/p].
- `report.json`: complete input snapshot, solver residual, FEM size, runtime,
  quadrature check and SCATMECH comparison when requested.
- `scatmech.npz`: reference calculated with lmax=14 and identical gap/materials.

FEM radiation amplitudes use the interface origin. Cross-solver complex phase
references have not been aligned; comparison currently uses absolute DSC only.
The angular hemisphere integral uses a 5 degree trapezoidal rule and periodic
azimuth average, equally for both solvers. Its own quadrature error is not yet
included in the reported discrepancy. All angular data are unnormalized DSC.

A small algebraic residual is NOT a physical error estimate. Check element
order/mesh, PML size/strength, domain size, angular integration, polarization
quadrature and contact-region field convergence separately. The two solvers
agreeing in far-field power does not prove the local contact field is accurate.
Near-field plots exclude an exact z=0 row; values very close to material
boundaries can still depend strongly on interpolation and mesh resolution.

Run analytic/native checks in the FEM environment:

```powershell
& '.\build\ngsolve_env\Scripts\python.exe' -m unittest discover -s tests -p test_fullwave.py -v
```

The full existing test suite runs in NIST2; the native FEM test is then skipped.

## Sources

- [NGSolve native installation](https://ngsolve.org/installation.html)
- [Netgen OCC materials and geometry](https://docu.ngsolve.org/ngs24/tutorials/03_occ.html)
- [NGSolve complex-coordinate PML](https://ngsolve.org/ngsolve/docs/i-tutorials/unit-1.7-helmholtz/pml.html)
- [NIST Bobbert-Vlieger reference](https://pages.nist.gov/SCATMECH/docs/bobvlieg.htm)

Next development gate: establish rough-patch mesh/domain/PML and near-field
convergence, independent roughness benchmarks and energy balance before
exposing this backend in the production UI. Ideal image-plane propagation is
available separately, but does not remove the FEM validation limitations.
