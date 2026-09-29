# tests/kernel/test_dynamic_parse_text_reporting_suite.m

## Purpose

Regression test suite for deterministic parsing, text, and safe reporting utilities in Spinach, covering operator-specification parsing, isotope predicates, label lookup, silent reporting calls, and polyadic text diagnostics.

## Behaviour

- Announces the test target with `fprintf` and initialises a regression result via `new_test_result` for `kernel/dynamic_parse_text_reporting_suite`.
- Builds a small three-spin system descriptor (`local_parse_spin_system`) with `sys.output` set to `'hush'`, isotopes `{'1H','E','13C'}`, labels `{'proton','electron','carbon'}`, all spin types `'S'`, multiplicities `[2 2 2]`, and coordinates `{[0 0 0],[1 0 0],[0 1 0]}`.
- Tests `human2opspec(spin_system,'Lz','nuclei')`: expects opspecs `{[2 0 0];[0 0 2]}` and coefficients `[1;1]`, verifying nuclei selection of non-electron spins and Lz mapping to IST index two.
- Tests `human2opspec(spin_system,{'Lx','Lz'},{1,3})`: expects product opspecs `{[1 0 2];[3 0 2]}` and Lx coefficients `[-sqrt(2);sqrt(2)]/2` (tolerances `1e-15`), reflecting the spherical tensor convention `Lx=(L+ + L-)/2`.
- Tests `idxof(sys,'carbon')` label lookup, expecting the one-based spin index `3`.
- Tests isotope predicates: `isnucleus('1H')` is true, `iselectron('E')` is true, and `isnucleus('E')` is false.
- Tests hush-mode reporting side-effect freedom using `evalc`: `report(spin_system,'hidden message')`, `banner(spin_system,'basis_banner')`, and `summary_coordinates(spin_system,'coordinate summary')` must all produce empty console text.
- Tests `polinfo(polyadic({{speye(2),sparse(1)}}))`: the captured text must contain `'polyadic [2x2]'`, `'matrix 1 [2x2]'`, and `'matrix 2 [1x1]'`, describing the unopened polyadic size and its matrix cores.

## Inputs and outputs

```matlab
result = test_dynamic_parse_text_reporting_suite()
```

- `result` — regression test result object with explanatory messages accumulated through `test_true` and `test_close` assertions.

## References

- [Source file on GitHub](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_dynamic_parse_text_reporting_suite.m)
