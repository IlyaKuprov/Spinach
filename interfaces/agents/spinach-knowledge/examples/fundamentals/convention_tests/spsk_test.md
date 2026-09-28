# examples/fundamentals/convention_tests/spsk_test.m

- Signature: `spsk_test()`

## Purpose

Checks the span-skew interaction convention implemented by `spsk2mat` against direct rotation of a diagonal tensor.

## Method and check

The test draws three separated eigenvalues (`xx=rand()`, `yy=rand()+3`, and `zz=rand()+6`) and random Euler angles. It constructs `AM=R*diag([xx yy zz])*R'` directly, then computes the isotropic value, span, and skew parameters used by `spsk2mat` to form `AS`. The 1-norm of `AM-AS` must be below 10⁻⁶.
