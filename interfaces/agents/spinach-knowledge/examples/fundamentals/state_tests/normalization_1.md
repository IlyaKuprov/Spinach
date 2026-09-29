# examples/fundamentals/state_tests/normalization_1.m

[Source code](https://github.com/IlyaKuprov/Spinach/blob/main/examples/fundamentals/state_tests/normalization_1.m)

## Purpose

This check asks whether the generated Cartesian angular-momentum states have equal norms for their x, y, and z components in each of Spinach's three listed formalisms. It checks component consistency within each representation; it does not compare the absolute state norms between formalisms.

## System and method

The test creates a two-spin system containing 1H and 235U at magnet field 14.1, with scalar Zeeman entries {2.5 1.0} and symmetric scalar-coupling entries of 10. For each of zeeman-hilb, zeeman-liouv, and sphten-liouv, it builds the basis with approximation='none', then constructs Lx, Ly, and Lz states for each isotope.

For each species and formalism, the script computes the three pairwise differences between MATLAB's default norm values after converting the states to full arrays. It collects six differences per formalism: three for 235U and three for 1H.

## Check and output

The check fails if any absolute norm difference exceeds 1e-6; otherwise it prints Internal norm consistency test PASSED. A failure raises Internal norm consistency test FAILED. The script reports only the pass/fail message, not the individual norms or residuals.

This is a small representation-consistency check for Cartesian state construction in the specified two-spin system. It does not test time evolution, other isotope combinations, reduced basis approximations, or physical observables.
