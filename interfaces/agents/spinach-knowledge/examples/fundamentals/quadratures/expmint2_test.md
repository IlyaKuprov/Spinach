# examples/fundamentals/quadratures/expmint2_test.m

## Purpose

This test compares the `expmint2` double-exponential matrix integral routine with a reference assembled from nested MATLAB numerical integrations. It checks one randomly drawn matrix problem per call; it is not a spin-system simulation.

## Test problem

The matrix dimension is `n = 10 + randi(5)`, hence 11 through 15, and the upper integration limit is `ul = randi(10) + rand()`, with `1 <= ul < 11`. Three independent random complex matrices are symmetrised to make A, C, and E Hermitian; B and D are general random complex matrices. The random stream is not seeded here, so different calls need not use the same case.

The example bootstraps a Spinach system and evaluates `expmint2(spin_system,A,B,C,D,E,ul)`. No other system or basis parameters are set in this file.

## Reference and decision rule

The reference nests two array-valued calls to MATLAB's `integral`. Its inner integrand is `exp(1i*C*x)*D*exp(-1i*E*x)`, with the inner result left-multiplied by `exp(-1i*C*t)`. The outer integrand is `exp(1i*A*t)*B*int_inner(t)`, and its integral is left-multiplied by `exp(-1i*A*T)`. Both integrations run from zero to their current upper variable, and the final reference is evaluated at `T = ul`.

The test forms the Frobenius-norm difference between the two outputs and divides it by the Spinach output's Frobenius norm. The comparison scale is written as the maximum of that same Spinach norm twice, so it reduces to that norm. The test raises an error only when this relative discrepancy is greater than `10*n*eps('double')`; otherwise it prints a pass message. That is the source's criterion, not a result reported here.

## Output and scope

A call prints either a pass message or raises the named failure error; it does not return a diagnostic value. The comparison is limited to a single unseeded random case, with MATLAB's default integration settings, and does not establish accuracy for other matrix structures, dimensions, upper limits, or composed applications. The source contains no reported run result or bibliographic citation.

Source: [expmint2_test.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/fundamentals/quadratures/expmint2_test.m).
