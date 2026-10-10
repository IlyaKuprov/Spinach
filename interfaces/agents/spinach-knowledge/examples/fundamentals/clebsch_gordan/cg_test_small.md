# examples/fundamentals/clebsch_gordan/cg_test_small.m

- MATLAB implementation: [examples/fundamentals/clebsch_gordan/cg_test_small.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/fundamentals/clebsch_gordan/cg_test_small.m)

- Signature: `cg_test_small()`

## Purpose

Checks Spinach Clebsch–Gordan coefficients against the arbitrary-precision Mathematica reference values recorded in the small test table. This is a table-driven numerical comparison, not a spectrum simulation.

## Test data and convention

Run `cg_test_small()` with `cg_test_table_small.mat` available in the MATLAB working directory and `clebsch_gordan` on the MATLAB path. The function loads the variable `cg_test_table`; for each row it passes columns 1–6, in order, to `clebsch_gordan` and treats column 7 as the reference value. The source does not give a separate interpretation or units for those six arguments, so their ordering is defined by the table and the call rather than by prose in this example.

## Comparison and output

For each row, the script times the coefficient call, computes the absolute difference from column 7, and prints the row and difference followed by `PASS` when the difference is strictly less than `2*eps`. A row at or above that threshold raises an error containing the difference and stops the function. The function has no return value; timing output is printed in seconds. No particular test outcome is asserted here.

The check covers only the cases present in `cg_test_table_small.mat`; the filename supplies no documented numerical size boundary or coverage guarantee. Agreement on this table alone does not establish accuracy for unlisted quantum-number combinations. The `small` distinction here is the selected table file, not a range specified by the MATLAB wrapper.
