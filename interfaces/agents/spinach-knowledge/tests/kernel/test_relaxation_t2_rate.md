# tests/kernel/test_relaxation_t2_rate.m

Source: [tests/kernel/test_relaxation_t2_rate.m](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_relaxation_t2_rate.m)

## Purpose

Regression test for the phenomenological T2 relaxation rate. It verifies that the `t1_t2` relaxation model assigns transverse `L+` order the negative generator eigenvalue corresponding to the specified R2 rate.

## Behaviour

- Announces the test target with `fprintf('TESTING: Phenomenological T2 decay rate\n')`.
- Initialises a regression test result via `new_test_result` under the identifier `kernel/relaxation_t2_rate`, with the description "Phenomenological T2 decay rate" and the specification "the t1_t2 model must assign transverse magnetisation the generator eigenvalue -R2.".
- Builds a one-spin system:
  - `sys.magnet = 14.1`
  - `sys.isotopes = {'1H'}`
  - `inter.zeeman.scalar = {0}`
  - `inter.relaxation = {'t1_t2'}`
  - `inter.r1_rates = {2.0}`
  - `inter.r2_rates = {7.0}`
  - `inter.equilibrium = 'zero'`
  - `inter.rlx_keep = 'secular'`
  - `inter.temperature = 298`
  - `bas.formalism = 'sphten-liouv'`
  - `bas.approximation = 'none'`
- Constructs the spin system with `test_spin_system(sys, inter, bas)`.
- Computes the relaxation superoperator `R = relaxation(spin_system)` and the transverse state `rho = state(spin_system, 'L+', '1H')`, then applies `Rrho = R * rho`.
- Checks that `L+` is an eigenstate with the negative R2 generator eigenvalue using `test_close(result, 'R2 eigenvalue on L+', Rrho, -7.0*rho, 1e-12, 1e-12, ...)`, with the message "decay generators have negative eigenvalues, with magnitude R2".

## Inputs and outputs

Syntax:

```matlab
result = test_relaxation_t2_rate()
```

- **Output**: `result` — regression test result structure with explanatory messages.
- **Input**: none.

## References

- [Spinach — tests/kernel/test_relaxation_t2_rate.m (GitHub)](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_relaxation_t2_rate.m)
