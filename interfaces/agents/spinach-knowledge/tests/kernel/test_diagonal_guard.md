# tests/kernel/test_diagonal_guard.m

## Purpose

Regression test for the diagonal relaxation retention formalism boundaries in Spinach. It verifies that unsupported Zeeman diagonal retention fails explicitly, while spherical-tensor diagonal retention and Zeeman full retention preserve trace and the specified relaxation rates. A shipped GISSMO XML fixture checks pure damping in both Liouville bases and the conversion from Lorentzian linewidth to transverse decay rate.

## Behaviour

- Builds spin-half (`1H`) and spin-one (`14N`) systems with a Lindblad bath (`inter.relaxation={'lindblad'}`, `inter.lind_r1_rates=4`, `inter.lind_r2_rates=7`, `inter.equilibrium='zero'`, `inter.temperature=298`, `inter.rlx_keep='labframe'`, `inter.rlx_dfs='keep'`, `sys.magnet=14.1`, `bas.approximation='none'`).
- For each isotope, in the `zeeman-liouv` formalism, builds the full Zeeman generator `R=relaxation(spin_system)` and checks:
  - Trace conservation: `unit'*R` equals `0*unit'` within `1e-12` (full retention preserves the left-null trace functional).
  - Identity stationarity: `R*unit` equals `0*unit` within `1e-12` (the symmetric Lindblad bath leaves identity stationary).
  - Omitted retention: removing the `keep` field from `spin_system.rlx` still yields the same full generator `R` (tolerances `0`).
  - Noncommuting fixture: with `bas_hilb.formalism='zeeman-hilb'`, a Hamiltonian `H=0.3*operator(spin_hilb,'Lx',1)+0.2*operator(spin_hilb,'Ly',1)` and a complex pure state `ket` with entries `1/sqrt(2)` and `1i/sqrt(2)`, forms `Q=hilb2liouv(H,'comm')` and requires `norm(Q*R-R*Q,'fro')>1e-3` (coherent and dissipative dynamics must not commute in this control).
  - Coherent evolution `rho_final=reshape(expm(full(-1i*Q+R))*rho(:),size(H))` must preserve trace (`trace(rho_final)` equals `1` within `1e-12`) and Hermiticity (`rho_final` equals `rho_final'` within `1e-12`).
  - Zeeman diagonal refusal: setting `spin_system.rlx.keep='diagonal'` and calling `relaxation` must throw an error whose message contains both `'not implemented'` and `'zeeman-liouv'`; both spin-half and spin-one diagonal Zeeman requests must be rejected explicitly.
- In the `sphten-liouv` formalism with `spin_system.rlx.keep='diagonal'`, checks:
  - Spherical trace conservation and identity stationarity within `1e-12`.
  - Longitudinal rate: `R*rho_z` equals `-4*rho_z` within `1e-12` (preserves the specified R1 rate `inter.lind_r1_rates=4`), where `rho_z=state(spin_system,'Lz',1)`.
  - Transverse rate: `R*rho_p` equals `-7*rho_p` within `1e-12` (preserves the specified R2 rate `inter.lind_r2_rates=7`), where `rho_p=state(spin_system,'L+',1)`.
- Damping-only case: with `inter.relaxation={'damp'}`, `inter.damp_rate=5`, and the Lindblad rate fields removed, for both `sphten-liouv` and `zeeman-liouv` the generator must equal `R_ref=-inter_damp.damp_rate*(unit_oper(spin_system)-unit*unit')` after `clean_up` with `spin_system.tols.rlx_zero`, within `1e-12` (non-selective decay with the unit state preserved).
- Imports a GISSMO subsystem from `examples/nmr_metabol/molecule_b.xml` with a one-hertz Lorentzian linewidth (second argument `2` to `gissmo2spinach`):
  - `inter.damp_rate` must equal `pi` exactly (tolerances `0`): the shipped one-hertz FWHM gives a pi-per-second damping rate.
  - `inter.rlx_keep` must be `'labframe'` (pure damping must not request unsupported Zeeman diagonal retention).
  - For both `sphten-liouv` and `zeeman-liouv`, the full generator `R` must match `R_ref=-pi*(unit_oper(spin_system)-unit*unit')` after `clean_up`, with `norm(R-R_ref,'fro')` equal to `0` within `1e-12`; the round bound is `spin_system.tols.rlx_zero*sqrt(nnz(R_ref))/2`.
  - GISSMO trace conservation (`unit'*R` vs `0*unit'`) and identity stationarity (`R*unit` vs `0*unit`) within `round_bound` and `1e-12` (entrywise rounding).
  - Transverse decay: with `rho_p=state(spin_system,'L+','1H')` and `rate=-real(rho_p'*R*rho_p)/(rho_p'*rho_p)`, requires `R*rho_p` equals `-pi*rho_p` within `round_bound*norm(rho_p)` and `1e-12` (the transverse signal decays as `exp(-pi*FWHM*t)`), and `rate/pi` equals `1` within `round_bound/pi` and `1e-12` (Lorentzian FWHM is the transverse decay rate divided by pi).
  - In `sphten-liouv` only, setting `spin_system.rlx.keep='diagonal'` and recomputing `R_diag=relaxation(spin_system)` must satisfy `norm(R-R_diag,'fro')` equal to `0` exactly (damping is added after retention, so the spherical result must not change).

## Inputs and outputs

```matlab
result = test_diagonal_guard()
```

- **Output**: `result` — regression result structure with explanatory messages, initialised via `new_test_result('kernel/diagonal_guard', 'Diagonal relaxation retention guard', 'unsupported Zeeman diagonal retention must fail without changing supported generators.')` and accumulated through `test_close` and `test_true` assertions.
- No inputs.

## References

- Source: [tests/kernel/test_diagonal_guard.m](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_diagonal_guard.m)
