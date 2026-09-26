# kernel/conventions/transforms/dcm2wigner.m

- Signature: `D=dcm2wigner(dcm)`

## Purpose

Converts a directional cosine matrix into second-rank Wigner function matrix. Syntax: D=dcm2wigner(dcm)

## Physical / mathematical content

- Convention and tensor-transform utilities. They convert among tensor parameterisations, coordinate systems, and unit systems; the underlying mathematics is linear algebra on rank-2 tensors and rotation representations.

## Numerical / algorithmic content

## Parameters / inputs

- dcm -directional cosine matrix

## Outputs

- D -matrix of second rank Wigner D functions. Rows
- and columns are sorted by descending ranks:
- [D( 2,2) ... D( 2,-2)
- ... ... ...
- D(-2,2) ... D(-2,-2)]
- Notes: the resulting Wigner matrix is to be used as v=W*v, where v is
- a column vector of irreducible spherical tensor coefficients in
- the following order: T(2,2), T(2,1), T(2,0), T(2,-1), T(2,-2).

## Implementation structure

- Converts a directional cosine matrix into second-rank Wigner function
- matrix. Syntax:
- D=dcm2wigner(dcm)
- dcm -directional cosine matrix
- D -matrix of second rank Wigner D functions. Rows
- and columns are sorted by descending ranks:
- [D( 2,2) ... D( 2,-2)
- ... ... ...
- D(-2,2) ... D(-2,-2)]
- a column vector of irreducible spherical tensor coefficients in
- the following order: T(2,2), T(2,1), T(2,0), T(2,-1), T(2,-2).
- Check consistency
