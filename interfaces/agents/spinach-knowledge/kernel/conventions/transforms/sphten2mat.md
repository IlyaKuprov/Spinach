# kernel/conventions/transforms/sphten2mat.m

- Signature: `M=sphten2mat(rank0,rank1,rank2)`

## Purpose

Converts the nine components of the irreducible spherical tensor re- presentation of an interaction tensor into the Cartesian representa- tion with a 3x3 matrix. The conventions are matched to Equation (18) of the paper by Len Mueller (http://dx.doi.org/10.1002/cmr.a.20224). Spherical tensor components should be listed in the following order: rank 0: (0,0) rank 1: (1,1) (1,0) (1,-1) rank 2: (2,2) (2,1) (2,0) (2,-1) (

## Physical / mathematical content

- Convention and tensor-transform utilities. They convert among tensor parameterisations, coordinate systems, and unit systems; the underlying mathematics is linear algebra on rank-2 tensors and rotation representations.

## Numerical / algorithmic content

## Parameters / inputs

- rank0 -a single number giving the coefficient of T(0,0) in
- the spherical tensor expansion.
- rank1 -a row vector with three numbers giving the coeffici-
- ents of T(1,1), T(1,0) and T(1,-1) in the spherical
- tensor expansion.
- rank2 -a row vector with five numbers giving the coeffici-
- ents of T(2,2), T(2,1), T(2,0), T(2,-1) and T(2,-2)
- in the spherical tensor expansion.

## Outputs

- M -3x3 interaction tensor

## Implementation structure

- Converts the nine components of the irreducible spherical tensor re-
- presentation of an interaction tensor into the Cartesian representa-
- tion with a 3x3 matrix. The conventions are matched to Equation (18)
- of the paper by Len Mueller (http://dx.doi.org/10.1002/cmr.a.20224).
- Spherical tensor components should be listed in the following order:
- rank 0: (0,0)
- rank 1: (1,1) (1,0) (1,-1)
- rank 2: (2,2) (2,1) (2,0) (2,-1) (2,-2)
- and should be supplied as coefficients in front of the corresponding
- irreducible spherical tensor operators returned by irr_sph_ten.m fun-
- ction Syntax:
- M=sphten2mat(rank0,rank1,rank2)
