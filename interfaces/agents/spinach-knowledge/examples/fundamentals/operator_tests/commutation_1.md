# examples/fundamentals/operator_tests/commutation_1.m

- Signature: `commutation_1()`

## Purpose

Checks spin-operator commutators across the three Spinach formalisms supported by this test: zeeman-hilb, zeeman-liouv, and sphten-liouv.

## System and checks

For each formalism, the script creates a single 1H spin at zero field with zero scalar chemical shift and uses basis approximation none. It constructs Lx, Ly, Lz, L+, and L−, then evaluates the Frobenius norms of the residuals for [Lz,L+]=L+, [Lz,L−]=−L−, and [Lx,Ly]=iLz.

The 3-by-3 array of residual norms passes if its Frobenius norm is below 1e-6. The script then reports that the cross-formalism test passed; otherwise it raises an error.
