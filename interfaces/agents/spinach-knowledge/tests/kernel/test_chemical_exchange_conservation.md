# tests/kernel/test_chemical_exchange_conservation.m

Source: [https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_chemical_exchange_conservation.m](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_chemical_exchange_conservation.m)

## Purpose

Regression test that verifies population conservation in a symmetric two-site chemical exchange model. The test checks that the kinetics generator built by `kinetics()` conserves the total population over the two sites, i.e. that probability leaving each site enters the other site.

## Behaviour

- Announces the test target with `fprintf('TESTING: Chemical exchange conservation\n')`.
- Registers a new test result via `new_test_result()` with the identifier `kernel/chemical_exchange_conservation`, the description `Chemical exchange conservation`, and the criterion string `closed two-site exchange must conserve total spin population.`.
- Builds a symmetric two-site exchange system with:
  - `sys.magnet = 14.1` (magnetic field, Tesla).
  - `sys.isotopes = {'1H','1H'}` (two protons).
  - `inter.zeeman.scalar = {0 0}`.
  - `inter.chem.parts = {1,2}` (two chemical exchange partitions).
  - `inter.chem.rates = [-3 3; 3 -3]` (symmetric rate matrix, s^-1).
  - `inter.chem.concs = [1 1]` (equal concentrations).
  - `bas.formalism = 'sphtten-liouv'` (spherical tensor Liouville space formalism).
  - `bas.approximation = 'none'`.
- Constructs the spin system with `test_spin_system(sys,inter,bas)`.
- Builds the kinetics generator with `K = kinetics(spin_system)`.
- Computes column sums of the full kinetics matrix: `col_sums = sum(full(K),1)`.
- Closes the test with `test_close(result,'zero column sums',col_sums,zeros(size(col_sums)),1e-15,1e-15,'probability leaving each site must enter the other site')`, asserting that all column sums equal zero within absolute and relative tolerances of `1e-15`.

## Inputs and outputs

```matlab
result = test_chemical_exchange_conservation()
```

- **Inputs**: none.
- **Outputs**:
  - `result` — regression test result structure with explanatory messages, as produced by `new_test_result()` and updated by `test_close()`.

## References

- Spinach `kinetics()` function — builds the chemical kinetics generator.
- Spinach `test_spin_system()` function — constructs a spin system from `sys`, `inter`, and `bas` specifications.
- Spinach testing framework functions `new_test_result()` and `test_close()`.
