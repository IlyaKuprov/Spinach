# tests/kernel/test_dynamic_remaining_core_suite.m

**Source:** https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_dynamic_remaining_core_suite.m

## Purpose

Regression test for the remaining deterministic utility helpers in the Spinach kernel. It verifies that small hand-written utility calls preserve documented algebraic and reporting semantics, covering block eliminations, text reporting, spin metadata, analytical line shapes, pumping terms, kite pruning, trajectory stitching, and small random-rotation diagnostics.

## Behaviour

The function announces its target with `fprintf`, registers a test result under the name `kernel/dynamic_remaining_core_suite` via `new_test_result`, and then runs a sequence of checks:

- **Adiabatic elimination** (`adelim`): builds a 3-spin Liouville-space descriptor with `local_liouvillian_system(3)` and applies `adelim` to the matrix `[1 2 3;4 5 6;7 8 10]` with state 3 eliminated and states `[1 2]` retained. The slow block must equal `[1 2;4 5]`, and the induced relaxation term must equal `1i*[3;6]*(1/10)*[7 8]`, i.e. `i*L01*(L11\L10)`, both at tolerances `1e-14`.
- **Clebsch-Gordan coefficient** (`cg_fast`): `cg_fast(1,0,1/2,1/2,1/2,-1/2)` must equal `1/sqrt(2)` at tolerance `1e-12` (two spin-half states coupling to the triplet M=0 state).
- **Console and banner reporting** (`report`, `banner`): writes to a temporary file handle stored in `spin_system.sys.output`, then verifies via `fileread` that the log contains both the direct report message and the text `BASIS SET`; `banner()` must delegate to `report()`. The output is then reset to `'hush'`.
- **Polyadic diagnostics** (`polinfo`): called directly on `polyadic({{speye(2),sparse(1)}})` with an explicit label; a PASS message is appended manually.
- **Coordinate summary** (`summary_coordinates`): called on a one-spin system with coordinates `[0 0 0]`, isotope `'1H'`, label `'proton'`, and multiplicity 2, traversing coordinate metadata through the silent report path.
- **Cell packaging** (`impound`): `impound(17,'spinach',{speye(2)})` must return a cell array whose entries equal the inputs unchanged.
- **Dipolar coupling** (`dipolar`): on a two-spin system with coordinates `[0 0 0]` and `[1 0 0]` (1 Angstrom apart on the X axis), gammas `[2 3]`, `hbar=1`, `mu0=4*pi`, the coupling tensor for pair `{1,2}` must equal `dip_pref*[-2 0 0;0 1 0;0 0 1]` where `dip_pref=0.5*gamma(1)*gamma(2)*hbar*mu0/(4*pi*(1e-10)^3)`; the proximity matrix must have exactly 2 nonzeros (both directed spin-pair proximities below the distance cutoff `prox_cutoff=2`).
- **Isotope swapping** (`isoswap`): on a two-spin system with isotopes `{'1H','13C'}`, bilinear couplings `eye(3)` and `2*eye(3)`, and a quadratic self-coupling `eye(3)` on spin 1, swapping spin 1 to `'2H'` must replace the isotope string, scale the `{1,2}` coupling by `spin('2H')/spin('1H')` to `gamma_ratio*eye(3)`, and wipe the quadratic self-coupling `{1,1}` (quadratic self-couplings are not transferable).
- **Interaction representation** (`intrep`): with `H0=2*pi*diag([0 1])` and `H=H0+diag([0.05 -0.02])`, order zero must return `H-H0=diag([0.05 -0.02])` after period validation, at tolerance `1e-13`.
- **Finite-difference Hessian** (`fdhess`): on a constant `ones(5,6,7)` field with a 3-point stencil, the output must be a `3x3` cell array and every block must equal `zeros(5,6,7)`; block inspection is gated on the layout guard passing.
- **Sinkhole removal** (`sinkhole`): on `L=[1 2 3;4 5 6;7 8 9]` with frozen states `[2 3]`, the result must be `[1 0 0;4 0 0;7 0 0]` with exact comparison (tolerances 0), zeroing all Liouvillian columns that feed frozen states.
- **Lorentzian line shape** (`lorentzcon`): the normalised branch with centre 0, amplitude `2*pi`, fwhm 2 on `x=[-1 0 1]` must equal `2./(1+x.^2)`; the narrow-linewidth branch with `fwhm=1e-308` at exact centre must equal `1/(pi*(fwhm/2))` without producing NaN from a reciprocal-width product; and, when a MEX version exists (`exist('lorentzcon','file')==3`), the segment branch at exact offsets `[0 1]` must return finite endpoint values `[0.5 0.5]` via two-argument `atan2`.
- **Gaussian line shape** (`gausscon`): the scalar-offset branch with centre 0, amplitude 2, width 2 must match `2*gaussfun(x,2)`; coincident triangle vertices `[0 0 0]` must collapse to the scalar-offset Gaussian; and a general triangle with offsets `[-1 0.25 2]`, amplitude 1, width 0.6 evaluated at `[-0.5 0.25 1.75]` must match direct numerical quadrature of the triangular weight convolved with `gaussfun`, using `integral` with waypoints and tolerances `1e-12`, compared at `1e-10`.
- **Magnetic pumping** (`magpump`): with `R=zeros(3)` and `rho=[0;2;-1]` at rate 0.25, the result must be `[0 0 0;0.5 0 0;-0.25 0 0]`, adding `rate*rho` to the first relaxation-superoperator column only.
- **Redfield-kite pruning** (`sec2kite`): on a 4-state system with basis `[0;1;2;3]` and a sparse relaxation matrix with entries at positions `(1,3)=2`, `(2,2)=7`, `(2,4)=5`, `(4,4)=9`, only the self-relaxation and longitudinal cross-relaxation entries `(1,3)=2`, `(2,2)=7`, `(4,4)=9` survive.
- **Sorensen bound** (`sorensen`): `sorensen(diag([3 1]),diag([2 0]))` must equal 1.5, the sorted-eigenvalue scalar product divided by `trace(rho_targ^2)`.
- **Trajectory stitching** (`stitch`): with zero Liouvillian, `rho_stack=[1 2;3 4]`, identity `coil_stack`, and step structures `t1.nsteps=2`, `t2.nsteps=3`, `t2.timestep=0.1`, `t3.nsteps=2`, every `t2` slice of the FID must equal `coil_stack'*rho_stack`.
- **Random rotations** (`rwalk`): with `rng(1,'twister')` and a one-worker process pool ensured by `local_ensure_pool`, `rwalk(5,1,1e-6)` must return a `5x3` array of finite Euler angles; when the shape guard holds, the first orientation must reconstruct the identity via `euler2dcm` (compared to `eye(3)` at `1e-14`).

## Inputs and outputs

```matlab
result=test_dynamic_remaining_core_suite()
```

**Inputs:** none.

**Outputs:**

- `result` — regression test result object with explanatory messages, created by `new_test_result` and accumulated through `test_close` and `test_true` checks.

## References

- Source file: [tests/kernel/test_dynamic_remaining_core_suite.m](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_dynamic_remaining_core_suite.m) in the Spinach repository.
- Functions exercised: `adelim`, `cg_fast`, `report`, `banner`, `polinfo`, `polyadic`, `summary_coordinates`, `impound`, `dipolar`, `isoswap`, `spin`, `intrep`, `fdhess`, `sinkhole`, `lorentzcon`, `gausscon`, `gaussfun`, `magpump`, `sec2kite`, `sorensen`, `stitch`, `rwalk`, `euler2dcm`, `new_test_result`, `test_close`, `test_true`.
