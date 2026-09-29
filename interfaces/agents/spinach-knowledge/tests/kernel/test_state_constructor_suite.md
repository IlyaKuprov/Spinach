# tests/kernel/test_state_constructor_suite.m

## Purpose

Regression test suite for Spinach state-constructor helper functions. It verifies that unit, thermal, singlet/triplet, partner, deuteron-pair, four-spin, and zero-field triplet state constructors produce correctly normalised physical density objects, using projector and normalisation identities.

Source: https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_state_constructor_suite.m

## Behaviour

- Announces the test target with `fprintf('TESTING: State-constructor functions\n')` and initialises a test result object via `new_test_result('kernel/state_constructor_suite', ...)`.
- **One-spin unit and thermal states**: builds a single `1H` spin system with zero scalar Zeeman interaction and `inter.temperature=300`.
  - `zeeman-hilb` formalism: `unit_state` must equal `speye(2)` (tolerances `1e-15`), and `equilibrium` with a zero Hamiltonian must equal `speye(2)/2` (tolerances `1e-14`).
  - `zeeman-liouv` formalism: `unit_state` must equal the stretched unit matrix `speye(2)` vectorised and divided by `sqrt(2)` (tolerances `1e-15`).
  - `sphten-liouv` formalism: `unit_state` must equal a sparse vector with 1 in the first position (the T(0,0) population basis vector, tolerances `1e-15`); also runs a `stateinfo(spin_s,unit_s,1)` smoke test.
- **Two-spin singlet/triplet projectors**: for two `1H` spins in `zeeman-hilb`, computes `S=singlet(spin2,1,2)` and `[TU,T0,TD]=triplet(spin2,1,2)`; checks unit traces of `S` and each triplet projector, completeness `S+TU+T0+TD=speye(4)`, singlet idempotence `S*S=S`, and orthogonality `trace(S*T0)=0` (all tolerances `1e-14`).
- **Partner states**: for three `1H` spins, calls `partner_state(spin3,{{'L+',2}},{{{'E','Lz'},[1 3]}})`; expects `numel(A)=4` combinations (exact match, tolerances 0) and verifies each returned state equals `state(spin3,descr{n},{1,2,3})` (tolerances `1e-14`).
- **Four-spin singlet-singlet state**: for four `1H` spins, checks `four_spin_states(spin4,1:4,'S(x)S')` equals `kron(singlet(spin2a,1,2),singlet(spin2a,1,2))` (tolerances `1e-14`).
- **Deuteron-pair states**: for two `2H` spins, calls `[Sd,Td,Qd]=deut_pair(spind,1,2)`; checks that the sum of singlet, triplet, and quintet projectors equals `speye(9)` and that all have unit traces (tolerances `1e-12`).
- **Zero-field triplet**: for an `E3` electron spin with `ZFS=diag([-1 0 1])*1e6` and `Z=eye(3)*28e9`, calls `zftrip(spine,ZFS,[0.2 0.3 0.5],Z,0.01,1)`; checks the resulting density matrix has unit trace and is Hermitian (tolerances `1e-12`).

## Inputs and outputs

```matlab
result = test_state_constructor_suite()
```

- **Outputs**: `result` — regression test result object with explanatory messages, accumulated through `test_close` and `test_true` comparisons.
- **Inputs**: none.

## References

- Source: https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_state_constructor_suite.m
