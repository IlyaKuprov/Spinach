# tests/kernel/test_dynamic_relaxation_models_suite.m

## Purpose

Regression test for the dynamic relaxation model helper functions in Spinach. The suite verifies that relaxation helpers assign mathematically expected rates on tiny spin systems, covering anisotropic and functional T1/T2 rates, damping, Lindblad relaxation, scalar Redfield, correlation functions, and the serial Redfield include.

## Behaviour

The test announces its target with `TESTING: Dynamic relaxation model helpers` and initialises a regression test result via `new_test_result` under the name `kernel/dynamic_relaxation_models_suite`.

It then performs the following checks, each appended to the result with explanatory messages:

- **Anisotropic T1/T2 tensor rates**: builds a one-spin `1H` system at 14.1 T with `inter.relaxation={'t1_t2'}`, diagonal R1 tensor `diag([1 3 5])` and R2 tensor `diag([2 4 6])`, `sphten-liouv` formalism with `approximation='none'`, and calls `rlx_t1_t2` at zero Euler angles. It verifies that the laboratory z axis samples the zz element of the R1 and R2 tensors, i.e. `r1_op*rho_z` equals `-5*rho_z` and `r2_op*rho_p` equals `-6*rho_p` (tolerances 1e-12).
- **Function-handle rates**: replaces the tensors with constant periodic function handles `@(alp,bet,gam)2+...` and `@(alp,bet,gam)3+...`, calls `rlx_t1_t2` at Euler angles `[0.2 0.3 0.4]`, and verifies the constant function handle sets the R1 and R2 rates directly (`-2*rho_z` and `-3*rho_p`, tolerances 1e-12).
- **Non-selective damping**: uses `inter.relaxation={'damp'}` with `damp_rate=5` and `rlx_keep='labframe'`, calls `relaxation`, and verifies the Liouville-space damping leaves the thermodynamic unit state unchanged while assigning rate `-5` to the `L+` state (tolerances 1e-12).
- **One-spin Lindblad rates**: uses `inter.relaxation={'lindblad'}` with `lind_r1_rates=4` and `lind_r2_rates=7`, and verifies the trace-preserving Lindblad generator leaves the unit state unchanged, damps `Lz` at rate `-4`, and damps `L+` at rate `-7` (tolerances 1e-12).
- **Scalar Redfield integral**: sets `spin_l.tols.rlx_integration=1e-5`, uses `H0=sparse(2,2)` (zero) and `H1=sparse([0 1; 1 0])` with `tau_c=1e-3`, calls `rlx_scalar` with `{[2 tau_c]}`, and compares against the analytic reference `-2*tau_c*(1-1e-5^2)*speye(2)` (tolerances 1e-11).
- **Correlation functions**: builds a Redfield-ready system with `sys_r.disable={'hygiene','asyredf'}`, `inter_r.relaxation={'redfield'}`, `tau_c={1e-9}`, `rlx_keep='labframe'`, `rlx_dfs='ignore'`. Calls `corrfun(spin_r,2,3,2,3,2)` and verifies the isotropic rank-two Wigner autocorrelation weight is `1/5`, the rate is `-1e9` (i.e. `-1/tau_c`), and the species state count is 3.
- **Axial rotational diffusion**: sets `spin_r.rlx.tau_c={[2e-9 5e-9]}` and verifies the diagonal weight remains `1/5` and the rate matches `-(6*D_eq+(2-2+1)^2*(D_ax-D_eq))` with `D_ax=1/(6*2e-9)` and `D_eq=1/(6*5e-9)` (rate tolerance 1e-6).
- **Rhombic rotational diffusion**: sets `spin_r.rlx.tau_c={[1e-9 2e-9 3e-9]}`, calls `corrfun(spin_r,2,3,3,3,3)`, and verifies the five analytical rank-two rates `-( [4*Dxx+Dyy+Dzz, Dxx+4*Dyy+Dzz, Dxx+Dyy+4*Dzz, 2*Dxx+2*Dyy+2*Dzz-2*delta, 2*Dxx+2*Dyy+2*Dzz+2*delta] )` where `Dxx=1/(6e-9)`, `Dyy=1/(12e-9)`, `Dzz=1/(18e-9)` and `delta=sqrt(Dxx^2+Dyy^2+Dzz^2-Dxx*Dyy-Dxx*Dzz-Dyy*Dzz)` (rate tolerance 1e-6).
- **Serial Redfield include**: sets `inter_r.zeeman.matrix={diag([1 2 -3])}` for an anisotropic one-spin system, calls `relaxation`, and verifies the resulting matrix is nonzero and preserves the thermodynamic unit state (tolerances 1e-12).

## Inputs and outputs

**Syntax**

```matlab
result = test_dynamic_relaxation_models_suite()
```

**Outputs**

- `result` - regression test result object with explanatory messages accumulated from each check.

The function takes no inputs.

## References

- Source: [tests/kernel/test_dynamic_relaxation_models_suite.m](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_dynamic_relaxation_models_suite.m)
