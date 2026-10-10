# examples/fundamentals/tensor_structures/kronm_test_1.m

- MATLAB implementation: [examples/fundamentals/tensor_structures/kronm_test_1.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/fundamentals/tensor_structures/kronm_test_1.m)

- Signature: `kronm_test_1()`
- Source: [`examples/fundamentals/tensor_structures/kronm_test_1.m`](../../../../../../examples/fundamentals/tensor_structures/kronm_test_1.m)

## Purpose

Checks that `kronm` applies an ordered list of Kronecker factors to a matrix of input columns consistently with explicitly forming the full Kronecker-product matrix and multiplying it by the same input. Both real and complex data are exercised.

## Mathematical scope and test data

This is a generic linear-algebra test, with no spin model, physical basis, or unit-bearing parameters. Each run chooses 3–6 square factors; every factor has side length 2–4. The input has 1–20 columns and `prod(dims)` rows. Real factors and input are random normal arrays; the complex test adds independent imaginary random components. The full reference `Q` is built by repeated Kronecker products of the factors.

## Use and checks

Run `kronm_test_1` in the Spinach MATLAB environment. For each data type the test compares `kronm(Q_terms,x)` with `Q*x`, accepting a one-norm difference below 1e-6; otherwise it raises an error. It prints the measured call durations for `kronm` and dense multiplication in seconds, but does not set a performance threshold. Those timing units are seconds; matrix sizes and tolerances are algebraic quantities without physical units.

The source's pass/error branches define test criteria only. They do not show that the test was executed or passed.
