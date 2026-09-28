# examples/fundamentals/tensor_structures/ttrain_test_1.m

- Signature: `ttrain_test_1()`

## Purpose

Checks `ttclass` arithmetic by evaluating `P*P+3*P` for a tensor train representing a three-factor Kronecker product, then comparing its materialised result with the same calculation on the explicit dense matrix.

## Physical / mathematical content

- The factors are `A=magic(5)`, a random `20x20` matrix `B`, and `C=1i*rand(15)`; each is normalised by its matrix 2-norm.
- `P=ttclass(1,{A;B;C},0)` represents `kron(A,kron(B,C))`. The tested expression is ordinary matrix multiplication and addition, not a spin-system calculation.

## Numerical / algorithmic content

- The tensor-train expression is converted with `full` and compared with the dense result using the one-norm criterion `norm(P_TT-P_US,1)<100*eps('double')`.

## Implementation structure

- Generate and normalise the three matrices, then create the tensor-train object.
- Evaluate the polynomial in `P` through `ttclass`, evaluate it again after forming the dense Kronecker product, and compare the results.
