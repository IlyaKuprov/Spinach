# tests/kernel/test_giant_ham_descr.m

**Source:** [https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_giant_ham_descr.m](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_giant_ham_descr.m)

## Purpose

Regression test for the giant spin Hamiltonian descriptor route. The test verifies that high-rank giant spin Hamiltonian terms assembled through the descriptor route match direct spherical-tensor assembly.

## Behaviour

- Announces the test target with `fprintf('TESTING: Giant spin Hamiltonian descriptor\n')` and initialises a regression test result via `new_test_result` under the identifier `kernel/giant_ham_descr`, with the message that high-rank giant spin terms must match direct spherical-tensor assembly.
- Runs two cases through the local helper `local_case`:
  - `'labframe'` assumption with `'strong'` strength and Euler angles `[0.41 0.29 0.13]`.
  - `'deer-zz'` assumption with `'secular'` strength and Euler angles `[0.17 0.39 0.51]`.
- Each case builds a compact high-rank giant spin system with `sys.magnet=0`, `sys.isotopes={'E8'}`, `inter.zeeman.scalar={0}`, giant-spin coefficient and Euler-angle cell arrays spanning ranks 1 through 3, `bas.formalism='zeeman-hilb'`, and `bas.approximation='none'`, using `test_spin_system` to construct the spin system.
- Applies the requested assumption with `assume`, builds the production Hamiltonian via `[I,Q]=hamiltonian(spin_system)` and orients it as `H_obs=I+orientation(Q,euler_angles)`, then builds a direct reference Hamiltonian `H_ref` via the local helper `local_ref`.
- Checks that high-rank components are present using `test_close` with the label `['giant spin rank depth ' strength]`, comparing `double(numel(Q)>=3)` against `1` with tolerances `0` and `0`, with the message that rank-three giant spin terms must extend the rotational basis.
- Checks the oriented Hamiltonian using `test_close` with the label `['giant spin Hamiltonian ' strength]`, comparing `H_obs` against `H_ref` with tolerances `1e-7` and `1e-12`, with the message that descriptor assembly must match direct spherical-tensor assembly.
- The reference builder `local_ref` starts from `mprealloc(spin_system,0)`, loops over spins and spherical ranks, computes Wigner matrices with `wigner(r,euler_angles(1),euler_angles(2),euler_angles(3))`, and handles the giant spin assumption:
  - `'strong'`: loops over spherical tensor projections `k=1:(2*r+1)`, contracts coefficients as `coeff=W(k,:)*spin_system.inter.giant.coeff{n}{r}(:)`, and adds `coeff*operator(spin_system,{ist_spec},{n})` with `ist_spec=['T' num2str(r) ',' num2str(r-k+1)]` when `abs(coeff)>spin_system.tols.liouv_zero`.
  - `'secular'`: contracts coefficients as `coeff=W(r+1,:)*spin_system.inter.giant.coeff{n}{r}(:)` and adds the `T<r>,0` operator contribution under the same tolerance check.
  - `'ignore'`: skips the term with `continue`.
- The reference Hamiltonian is Hermitised as `H_ref=(H_ref+H_ref')/2` to match `orientation.m`.

## Inputs and outputs

```matlab
result = test_giant_ham_descr()
```

- **Output:** `result` — regression test result with explanatory messages.
- No inputs are required.

## References

- Source file: [tests/kernel/test_giant_ham_descr.m](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_giant_ham_descr.m) in the Spinach repository.
