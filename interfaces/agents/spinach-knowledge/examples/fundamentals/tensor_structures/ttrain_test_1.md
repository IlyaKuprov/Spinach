# examples/fundamentals/tensor_structures/ttrain_test_1.m

- Signature: `ttrain_test_1()`
- Source: [`examples/fundamentals/tensor_structures/ttrain_test_1.m`](https://github.com/IlyaKuprov/Spinach/blob/main/examples/fundamentals/tensor_structures/ttrain_test_1.m)

## Purpose

A no-argument MATLAB test of arithmetic on a `ttclass` tensor-train object. It checks the materialised result of `P*P+3*P` against the same polynomial evaluated on the explicit dense Kronecker-product matrix.

## Tensor factors and computation

The factors are `A=magic(5)`, `B=randn(20)`, and `C=1i*rand(15)`; each is divided by its matrix 2-norm. The source constructs `P=ttclass(1,{A;B;C},0)`, then computes `P_TT=full(P*P+3*P)`. The dense reference is formed in the same factor order, `kron(A,kron(B,C))`, and the expression `P_US=P*P+3*P` is evaluated on that matrix. The factor dimensions give a 1500×1500 reference matrix.

This is a generic matrix-arithmetic check, not a physical spin model: the source specifies no spins, Hamiltonian, basis, or physical units. The constructor arguments are recorded as written; this test does not explain their meanings. The random factors are unseeded.

## Comparison, output, and limits

The check is the strict one-norm condition `norm(P_TT-P_US,1)<100*eps('double')`; `eps('double')` is MATLAB's double-precision machine epsilon. If the condition is true, the source displays `Test passed.`; otherwise it raises `Test failed.`. The function declares no output arguments. This compares one polynomial on one randomly generated three-factor example; it is not evidence of an executed pass, a performance measurement, or validation for other tensor trains.
