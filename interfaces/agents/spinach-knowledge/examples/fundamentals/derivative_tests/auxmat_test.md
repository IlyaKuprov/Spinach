# examples/fundamentals/derivative_tests/auxmat_test.m

- MATLAB implementation: [examples/fundamentals/derivative_tests/auxmat_test.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/fundamentals/derivative_tests/auxmat_test.m)

- Signature: `auxmat_test()`

## Purpose

Checks the auxiliary block-matrix identity for differentiating the matrix exponential. For a matrix-valued `L(a,b,c)` and one parameter direction `E=dL/dt`, the upper-right `dim×dim` block of

`expm([L E; zeros(dim) L])`

is compared with a finite-difference derivative of `expm(L)`. In the source, `f=@(A)expm(A)`; the test is about the exponential, not an arbitrary matrix function.

## Test inputs and parameterisation

A random dimension `dim=randi(20)` is used. The complex random matrices `A`, `B`, and `C` are normalised in the spectral norm; as written, the normalisation of `C` divides by `norm(B,2)`. The matrix family is

`L(a,b,c)=a*A+c*cos(b)*B+(a*c^4)*C+a*b*c*B*C`.

The source uses the partial derivatives `dL/da=A+c^4*C+b*c*B*C`, `dL/db=-c*sin(b)*B+a*c*B*C`, and `dL/dc=cos(b)*B+4*a*c^3*C+a*b*B*C`. Independent real random values `x`, `y`, and `z` are used for the three checks, with `h=sqrt(eps)`.

## Comparison and acceptance criterion

For each parameter, the auxiliary block `P` is obtained from the block matrix exponential using that parameter's derivative of `L`. The reference `Q` is a fourth-order centred finite difference of `expm(L)`:

`Q = 2*(f(t+h)-f(t-h))/(3*h) - (f(t+2*h)-f(t-2*h))/(12*h)`,

with the other two parameters held fixed. The test raises `differentiation test failed` when `norm(P-Q,2)/norm(P,2)>1e-3`; otherwise it prints that parameter's pass message. The comparison is repeated separately for all three parameters.

## Scope

The dimension, matrices, and evaluation point are random, and the source sets no random seed. The script checks one point per parameter. Its relative-error denominator is `norm(P,2)` without an explicit zero-denominator guard.
