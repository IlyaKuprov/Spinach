# kernel/conventions/transforms/mat2sphten.m

Source implementation: https://github.com/IlyaKuprov/Spinach/blob/main/kernel/conventions/transforms/mat2sphten.m
Spinach Wiki: [mat2sphten.m](https://spindynamics.org/wiki/index.php?title=mat2sphten.m)
Reference: Len Mueller, Equation (22), [DOI: 10.1002/cmr.a.20224](https://doi.org/10.1002/cmr.a.20224).

## Purpose and convention

Converts a real 3x3 interaction matrix into coefficients of irreducible spherical tensor operators: one rank-0, three rank-1, and five rank-2 components. The convention follows Equation (22) of Len Mueller. Components are ordered by decreasing q within each rank: rank 1 as T(1,1), T(1,0), T(1,-1), and rank 2 as T(2,2), T(2,1), T(2,0), T(2,-1), T(2,-2). No coordinate rotation is applied; the matrix indices are consumed in the supplied ordering. The outputs are coefficients for the irreducible spherical-tensor operators returned by [irr_sph_ten.m](../../operators/irr_sph_ten.md).

## Syntax

```matlab
[rank0,rank1,rank2]=mat2sphten(M)
```

## Input

- `M` is a real numeric 3x3 interaction matrix. An empty input is replaced by `zeros(3)` before validation, so it produces zero coefficients. Nonempty inputs must be exactly 3x3.

## Outputs

- `rank0` is the scalar coefficient of T(0,0): `trace(M)/3`.
- `rank1` is a 3x1 column vector containing the coefficients of T(1,1), T(1,0), and T(1,-1), in that order.
- `rank2` is a 5x1 column vector containing the coefficients of T(2,2), T(2,1), T(2,0), T(2,-1), and T(2,-2), in that order.

For input elements M(i,j), the components are computed as:

- `rank1(1)=-(1/2)*(M(3,1)-M(1,3)-1i*(M(3,2)-M(2,3)))`
- `rank1(2)=+(1i/sqrt(2))*(M(1,2)-M(2,1))`
- `rank1(3)=-(1/2)*(M(3,1)-M(1,3)+1i*(M(3,2)-M(2,3)))`
- `rank2(1)=+(1/2)*(M(1,1)-M(2,2)-1i*(M(1,2)+M(2,1)))`
- `rank2(2)=-(1/2)*(M(1,3)+M(3,1)-1i*(M(2,3)+M(3,2)))`
- `rank2(3)=+(1/sqrt(6))*(2*M(3,3)-M(1,1)-M(2,2))`
- `rank2(4)=+(1/2)*(M(1,3)+M(3,1)+1i*(M(2,3)+M(3,2)))`
- `rank2(5)=+(1/2)*(M(1,1)-M(2,2)+1i*(M(1,2)+M(2,1)))`

The components are linear combinations of the supplied matrix entries; the output coefficients may be complex even though `M` must be real.
