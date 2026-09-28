# kernel/overloads/@ttclass/mldivide.m

- Signature: `x=mldivide(A,y)`

## Purpose

Computes a tensor-train solution to a linear system using the AMEn solver.

## Numerical / algorithmic content

The implementation shrinks the operands, forms the symmetrised system `(A'*A)*x=A'*y), and calls `amensolve` with tolerance `1e-6`. This is a least-squares normal-equation solve, not a direct unsymmetrised solve.

## Parameters / inputs

- A -ttclass matrix
- y -ttclass vector

## Outputs

- x -ttclass vector
