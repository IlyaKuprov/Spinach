# kernel/overloads/@ttclass/sum.m

- Signature: `answer=sum(ttrain,dim)`

## Purpose

Sum of elements of a tensor train representation of a matrix. Syntax: answer=sum(ttrain,dim)

## Physical / mathematical content

- Tensor-train linear algebra. These files implement compressed high-dimensional operators and AMEn/SVD-based algebra in tensor-train format.
- The relevant state manifold is the singlet/triplet decomposition, where permutation symmetry controls selection rules, relaxation susceptibility, and convertibility to ordinary magnetisation.

## Numerical / algorithmic content

## Parameters / inputs

- ttrain -tensor train object representing
- a matrix
- dim -summation dimension, 1 or 2

## Outputs

- answer -tensor train or flat representation
- of the summation result

## Implementation structure

- Sum of elements of a tensor train representation of
- a matrix. Syntax:
- answer=sum(ttrain,dim)
- ttrain -tensor train object representing
- a matrix
- dim -summation dimension, 1 or 2
- answer -tensor train or flat representation
- of the summation result
- Get sizes and ranks
- If all dimensions are singleton, return a scalar immediately
- In dim is omitted, choose first non-singleton dimension
- (this mimics the Matlab behaviour for matices)
