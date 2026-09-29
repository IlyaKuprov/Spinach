# examples/fundamentals/state_tests/normalization_2.m

[Source code](https://github.com/IlyaKuprov/Spinach/blob/main/examples/fundamentals/state_tests/normalization_2.m)

## Purpose

This check compares state inner products across Spinach's Zeeman Hilbert, Zeeman Liouville, and spherical-tensor Liouville formalisms, including the normalisation factor used by the source for the spherical-tensor representation. Unlike normalisation_1, it tests squared norms rather than differences between Cartesian-component norms.

## System and method

The script uses a two-spin system with 1H and 235U, magnet field 14.1, scalar Zeeman entries {2.5 1.0}, and symmetric scalar-coupling entries of 10. For each of zeeman-hilb, zeeman-liouv, and sphten-liouv, it uses approximation='none' and constructs Lx, Ly, and Lz states for both spins.

For each of the six states, it evaluates the squared Frobenius norm as `trace(A'*A)`, where A is the full state array. These values form six rows, in x/y/z order for 235U followed by x/y/z for 1H, and three columns in the formalism order above. The scaling factor is the product of the spin multiplicities in `spin_system.comp.mults`.

## Check and output

The code compares the Zeeman-Hilbert and Zeeman-Liouville columns directly, and compares the Hilbert column with the spherical-tensor column multiplied by the multiplicity product. Each comparison uses the L1 norm of the six-element difference and fails when it exceeds 1e-6. On failure it raises Cross-formalism state norm test FAILED; otherwise it prints Cross-formalism state norm test PASSED. The internally assembled six-by-three array is not printed by the script.

The conclusion is limited to these six constructed states, this two-spin configuration, and the three untruncated formalisms. This example does not test dynamics or establish equality without the spherical-tensor multiplicity scaling used in its comparison.
