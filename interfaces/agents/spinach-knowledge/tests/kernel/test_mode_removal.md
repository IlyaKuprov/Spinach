# tests/kernel/test_mode_removal.m

**Source:** [https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_mode_removal.m](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_mode_removal.m)

## Purpose

Regression test for `kill_spin` particle removal. It verifies that removing spins or bosonic modes from a Spinach system preserves the Hamiltonians and relaxation superoperators of the retained particles, as checked against independently constructed reference systems.

## Behaviour

- Creates a test result object via `new_test_result('kernel/mode_removal', ...)` describing retained mode reindexing.
- Iterates over ten hit lists covering spectator, spin, mode, simultaneous, logical-index, and no-op removals: `2`, `3`, `1`, `[2 4]`, `[false true false false false]`, `5`, `[]`, `[1 5]`, `[2 3 4]`, `[1 2 3 4]`.
- For each case, builds a five-particle system (`C3`, `13C`, `2H`, `1H`, `V3`) with two bosonic modes, removes the specified particles with `kill_spin`, and compares against a freshly built system containing only the retained particles.
- Checks that `trimmed.inter.modes` equals the reference mode data exactly, including nested spin leaves; when all modes are removed, checks that no `modes` container remains.
- For both `zeeman-hilb` and `zeeman-liouv` formalisms, applies the `labframe` assumption, adds an orientation term with direction `[0.17 0.31 0.23]`, and compares Hamiltonians with tolerances `1e-10` (absolute) and `1e-12` (relative).
- For the first case in the first formalism, asserts that the reference Hamiltonian has Frobenius norm of the imaginary part greater than 1 (transverse y fields produce complex matrix elements), that the Hamiltonian does not commute with the mode number operator `operator(expected,'N',1)` (noncommuting mode terms), and that the Hamiltonian is Hermitian.
- In the Liouville formalism, compares `relaxation` superoperators between trimmed and reference systems with the same tolerances, and asserts a nonzero reference decay norm when modes remain.
- For bosonic particle types `C3`, `V3`, and `T3`, builds mixed spin-mode systems with `1H` and `13C` spins (Zeeman scalars 1 and 2, mode frequency and carrier 1000) and tests assumption resets: after `kill_spin` removes the mode, `trimmed.inter.modes` must no longer contain a `strength` field.
- Rebuilds each system under `labframe`, `cavity`, and `spin-phonon` assumptions in both formalisms and compares Hamiltonians against an independently built two-particle spin-mode reference.
- Removes the final mode from four sources (unassumed, `labframe`, `cavity`, `spin-phonon` assumed), verifies the `modes` container is gone, and compares spin-only Hamiltonians against an independent `1H`/`13C` reference under `nmr` and `esr` assumptions with retention options `''`, `'zeeman'`, and `'couplings'`.
- The `build_system` helper constructs systems directly with `create` rather than using the production removal path, assigning mode frequencies 800 and 1300, carriers 10000 and 17000, anharmonicities 13 and 29, linewidths 11 and 19, and `t2_times` 0.003 and 0.007 to the two mode particles, plus pair couplings (`exchange`, `kerr`, `longitudinal`, `dispersive`) with strengths 17, 19, 23, and 29, spin-one quadratic tensors, Zeeman vectors, and first-derivative (`coupling_mod`, `zeeman_mod`) Raman terms.
- The `grumble` helper enforces that the keep list is an increasing row vector of particle identities from 1 to 5.

## Inputs and outputs

**Syntax:**

```matlab
result = test_mode_removal()
```

**Outputs:**

- `result` — test result object accumulating regression checks against independently created systems.

**Inputs:**

- None.

## References

- Uses `new_test_result`, `kill_spin`, `assume`, `basis`, `hamiltonian`, `orientation`, `operator`, `relaxation`, `create`, `test_true`, and `test_close`.
