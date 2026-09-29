# kernel/utilities/tolerances.m

Source: [kernel/utilities/tolerances.m](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/utilities/tolerances.m)

## Purpose

Sets the accuracy cut-offs, constants and tolerances used by the Spinach kernel. The function parses the `sys.tols` substructure into `spin_system.tols`, applies safe defaults, paranoid or loose presets depending on the `sys.enable` switches, and writes fundamental physical constants into the system description object.

## Behaviour

- Syntax: `[spin_system,sys]=tolerances(spin_system,sys)`.
- The function first calls its internal consistency checker `grumble(spin_system,sys)`, which validates every user-specified `sys.tols` field (numeric type, realness, scalar nature, sign, integrality or range as appropriate) and errors with a descriptive message on violation.
- For each tolerance, a user-specified value in `sys.tols` takes precedence; the field is then removed from `sys.tols`. Otherwise the value is chosen from the `paranoia` or `cowboy` presets when the corresponding string is present in `spin_system.sys.enable`, and otherwise from the safe default. Each setting is reported via `report(spin_system,...)` with the tag `user-specified`, `paranoid`, `loose` or `safe default`.
- Interaction tensor clean-up tolerance `inter_cutoff` (Hz): `eps()` when paranoid, `1e-2` when cowboy, `1e-10` safe default.
- Liouvillian matrix element zero tolerance `liouv_zero` (rad/s): `eps()` when paranoid, `1e-5` when cowboy, `1e-10` safe default.
- Relaxation superoperator element zero tolerance `rlx_zero` (Hz): `eps()` when paranoid, `1e-5` when cowboy, `1e-10` safe default.
- Propagator matrix element zero tolerance `prop_chop`: `eps()` when paranoid, `1e-8` when cowboy, `1e-10` safe default.
- Subspace population tolerance `subs_drop`: `eps()` when paranoid, `1e-2` when cowboy, `1e-10` safe default.
- Irrep population tolerance `irrep_drop`: `eps()` when paranoid, `1e-2` when cowboy, `1e-10` safe default.
- Steady state convergence tolerance `stst_tol`: `eps()` when paranoid, `1e-6` when cowboy, `1e-8` safe default.
- ZTE sample length `zte_nsteps`: `16` when cowboy, `NaN` (ZTE disabled) when paranoid, `32` safe default.
- ZTE zero track tolerance `zte_tol`: `1e-6` when cowboy, `NaN` (ZTE disabled) when paranoid, `1e-24` safe default.
- ZTE state vector density threshold `zte_maxden`: `NaN` (ZTE disabled) when paranoid, `0.5` safe default; there is no cowboy branch.
- Proximity tolerance for dipolar couplings `prox_cutoff` (Angstrom): `inf()` when paranoid, `3.5` when cowboy, `100` safe default. The source carries a TODO note to replace this with an energy tolerance.
- Krylov method switchover `krylov_tol`: `10000` safe default; no paranoid or cowboy branches.
- Basis printing hush tolerance `basis_hush`: `256` safe default; no paranoid or cowboy branches.
- Subspace bundle size `merge_dim`: `1000` safe default; no paranoid or cowboy branches.
- Sparse algebra tolerance on density `dense_matrix`: `0.15` safe default; no paranoid or cowboy branches.
- Sparse algebra tolerance on dimension `small_matrix`: `200` safe default; no paranoid or cowboy branches.
- Relative accuracy of the elements of the Redfield superoperator `rlx_integration`: `1e-6` when paranoid, `1e-2` when cowboy, `1e-4` safe default.
- Algorithm selection for propagator derivatives `dP_method`: default string `'auxmat'`; no paranoid or cowboy branches.
- Number of PBC images for dipolar couplings `dd_ncells`: `2` safe default; no paranoid or cowboy branches.
- Cache storage timeout `cache_mem` (days before a cache record is deleted): `365` safe default; no paranoid or cowboy branches.
- Fundamental constants are then written unconditionally: `hbar = 6.62607015e-34/(2*pi)` (J*s, exact number), `kbol = 1.380649e-23` (J*K^-1, exact number), `freeg = 2.00231930436092` (CODATA 2022), `mu0 = 1.25663706127e-6` (N A^-2, CODATA 2022) and `muB = 9.2740100657e-24` (J T^-1, CODATA 2022).
- If `paranoia` is enabled, the function appends `'zte'` to `spin_system.sys.disable` (so zero track elimination is disabled) and removes `'op_cache'`, `'prop_cache'` and `'ham_cache'` from `spin_system.sys.enable` (so operator, Hamiltonian and propagator caching are not enabled).
- After parsing, any fields remaining in `sys.tols` are reported as `unrecognised option` and the function errors with `unrecognised options in sys.tols`; the `tols` field is then removed from `sys`.
- The `grumble` consistency checker also validates a `path_drop` field that the main function body does not parse; it must be a non-negative real scalar if supplied.
- The header notes that direct calls and modifications to this function are discouraged: accuracy settings should be modified by setting the `sys.tols` structure, as described in the input preparation manual.

## Inputs and outputs

**Inputs**

- `spin_system` — Spinach system description object.
- `sys` — system specification object described in the input preparation section of the manual; may carry a `tols` substructure with user-specified tolerance fields.

**Outputs**

- `spin_system` — updated system description object with the populated `spin_system.tols` structure.
- `sys` — system specification structure with the tolerance substructure parsed out.

## References

- [Spinach Wiki: tolerances.m](https://spindynamics.org/wiki/index.php?title=tolerances.m)
- [Source file on GitHub](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/utilities/tolerances.m)
