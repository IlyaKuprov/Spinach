# tests/kernel/test_dynamic_physics_formulae_suite.m

**Source:** [https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_dynamic_physics_formulae_suite.m](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_dynamic_physics_formulae_suite.m)

## Purpose

Regression test for the deterministic closed-form physical formula utility helpers in the Spinach kernel. The suite verifies that spin addition, point-dipole tensors, hyperfine tensors, exponential drops, skew-normal densities, oscillator grids, hydrodynamic derivative construction, and spherical-tensor projection metadata match their analytic definitions on small systems.

## Behaviour

The function announces the test target with `fprintf('TESTING: Physical formula utilities\n')` and initialises a test result object via `new_test_result` for the target `kernel/dynamic_physics_formulae_suite`, described as "Physical formula utilities" with the property that closed-form physical helper formulae must match their analytic definitions on small systems.

The suite then performs the following checks:

- **Spin addition (`add_spins`)**: Called with two spin-half irreps (`add_spins(1/2,1/2)`). Verifies that the returned multiplicities are `[1 3]` (one singlet and one triplet irrep). With the projectors concatenated as `P=[proj{1} proj{2}]`, checks that `P'*P` equals `eye(4)` and `P*P'` equals `eye(4)`, both to tolerances `1e-12`, confirming the canonical projectors are orthonormal and span the four-dimensional product space.
- **Point-dipole coupling (`xyz2dd`)**: Called for a one-Angstrom z-axis displacement between `1H` and `13C` at origins `[0 0 0]` and `[0 0 1]`. Using `hbar=1.054571628e-34` and `mu0=4*pi*1e-7`, the reference coupling is `d_ref=spin('1H')*spin('13C')*hbar*mu0/(4*pi*(1e-10)^3)`. Checks that the coupling constant `d` matches `d_ref` (tolerances `1e-6`, `1e-12`), that the Euler angles `[alp bet gam]` are `[0 0 0]` (tolerances `1e-14`), and that the tensor `M` equals `d_ref*diag([1 1 -2])` (tolerances `1e-6`, `1e-12`), i.e. traceless principal values `d`, `d`, `-2d`.
- **Point hyperfine tensor (`xyz2hfc`)**: Called for a z-axis displacement with `1H`. With `C=1e4*spin('1H')*hbar*mu0/(4*pi*(1e-10)^3)`, checks that the returned tensor `A` equals `C*diag([-1 -1 2])` to tolerances `1e-4`, `1e-12` (a Gauss hyperfine tensor proportional to `diag(-1,-1,2)`).
- **Exponential drop (`expdrop`)**: Called as `expdrop(5,1,2,3,rate)` with `rate=log(4)`. Checks that the result equals `[5 9/5 1]` to tolerances `1e-14`, i.e. with rate `log(4)` over two seconds the middle exponential factor is one quarter.
- **Skew-normal density (`snormpdf`)**: Evaluated at `x=[-1 0 1]` with zero mean, unit standard deviation, and zero skew. Checks that the result equals the ordinary normal density `exp(-(x.^2)/2)/sqrt(2*pi)` to tolerances `1e-15` (Azzalini skew-normal reduces to the normal density when alpha is zero).
- **Oscillator grid and operators (`oscillator`)**: Called with `parameters.frc_cnst=4`, `parameters.par_mass=2`, `parameters.grv_cnst=0`, `parameters.n_points=11`, `parameters.box_size=10`. Checks that the grid `xgrid` equals `(-5:5)'` (tolerances `1e-15`), that `diag(X_oscl)` equals `xgrid` (tolerances `1e-15`), and that `H_oscl-H_oscl'` equals `zeros(11)` (tolerances `1e-15`), confirming the finite-difference Hamiltonian is Hermitian.
- **Hydrodynamic derivatives (`hydrodynamics`)**: Called with `hydro_spin_system.sys.enable={'polyadic'}`, `hydro_params.dims=10`, `hydro_params.npts=10`, `hydro_params.deriv={'period',3}`. With `Dx=fdmat(10,3,1)/(hydro_params.dims/hydro_params.npts)`, checks that `inflate(Fx)` equals `-1i*Dx` to tolerances `1e-15` (the one-dimensional derivative wrapped in a polyadic), and that `Fy` and `Fz` are both empty.
- **Spherical-tensor to Zeeman projection (`sphten2zeeman`)**: Called with `spin_system.comp.mults=2`, `spin_system.bas.formalism='sphten-liouv'`, `spin_system.bas.basis=(0:3)'`. Checks that `size(P)` is `[4 4]` and `rank(full(P))==4`, i.e. one spin-half has four spherical-tensor basis states and four Zeeman-Liouville states.

## Inputs and outputs

```matlab
result = test_dynamic_physics_formulae_suite()
```

**Outputs**

- `result` — regression test result object with explanatory messages, accumulated through the `test_true` and `test_close` assertions described above.

**Inputs**

- None.

## References

- Source file: [tests/kernel/test_dynamic_physics_formulae_suite.m](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_dynamic_physics_formulae_suite.m) in the Spinach repository.
- Functions exercised: `add_spins`, `xyz2dd`, `xyz2hfc`, `expdrop`, `snormpdf`, `oscillator`, `hydrodynamics`, `sphten2zeeman`, with supporting utilities `new_test_result`, `test_true`, `test_close`, `spin`, `fdmat`, `inflate`.
