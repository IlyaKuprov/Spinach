# tests/kernel/test_orientation_average_suite.m

Source: [tests/kernel/test_orientation_average_suite.m](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_orientation_average_suite.m)

## Purpose

Regression test suite for the `orientation()` and `average()` helper functions, exercising them on small exact cases. It checks the zero-Euler-angle orientation contraction against an explicit diagonal Wigner sum, and verifies that first-order average Hamiltonian theory leaves an unmodulated Hamiltonian unchanged.

## Behaviour

- Announces the test target with `fprintf('TESTING: Orientation and average helpers\n')`.
- Initialises a regression test result via `new_test_result('kernel/orientation_average_suite', ...)`, describing the target as rotational helper kernels that must preserve exact limiting cases.
- Builds a synthetic rank-one rotational basis `Q` as a 3-by-3 cell array of 2-by-2 sparse matrices, with diagonal components `Q{1}{1,1} = [1 0; 0 0]`, `Q{1}{2,2} = [0 1; 1 0]`, and `Q{1}{3,3} = [0 0; 0 2]`; all off-diagonal components are zero sparse matrices.
- Contracts the zero-orientation Hamiltonian as `H = orientation(Q,[0 0 0])` and compares it to the reference sum `H_ref = Q{1}{1,1} + Q{1}{2,2} + Q{1}{3,3}` using `test_close` with tolerances `1e-14` (absolute and relative), under the name `zero Euler orientation`. The rationale recorded is that `wigner(r,0,0,0)` is the identity, so only diagonal rotational components remain.
- Builds a quiet spin system for `average()` diagnostics: `sys.magnet = 0`, `sys.isotopes = {'1H'}`, `inter.zeeman.scalar = {0}`, `bas.formalism = 'sphten-liouv'`, `bas.approximation = 'none'`, passed through `test_spin_system`.
- Defines an unmodulated Hamiltonian decomposition with `Hp` and `Hm` as zero 2-by-2 sparse matrices, `H0 = [0 1; -1 0]` sparse, and modulation frequency `omega = 2*pi*1000`.
- Runs `H_avg = average(spin_system,Hp,H0,Hm,omega,'ah_first_order')` and compares it to `H0` using `test_close` with tolerances `1e-14` (absolute and relative), under the name `unmodulated average`. The rationale recorded is that with zero positive and negative Fourier components, first-order averaging returns `H0`.

## Inputs and outputs

Syntax:

```matlab
result = test_orientation_average_suite()
```

- `result` — regression test result structure with explanatory messages, accumulated through `test_close` checks.
- The function takes no input arguments.

## References

- `orientation` — rotational contraction helper under test.
- `average` — average Hamiltonian theory helper under test.
- `new_test_result`, `test_close`, `test_spin_system` — test harness utilities used by the suite.
