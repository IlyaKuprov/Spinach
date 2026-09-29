# kernel/conventions/transforms/sphten2mat.m

[MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/conventions/transforms/sphten2mat.m) · [Spin Dynamics Wiki](https://spindynamics.org/wiki/index.php?title=sphten2mat.m) · [Len Mueller paper, Equation (18), DOI](https://doi.org/10.1002/cmr.a.20224)

- Signature: `M=sphten2mat(rank0,rank1,rank2)`

## Purpose and component order

Convert coefficients of the irreducible spherical tensor operators returned by `irr_sph_ten.m` into a `3x3` Cartesian interaction tensor, using the convention matched to Equation (18) of Len Mueller’s paper. Supply components in this exact order:

- `rank0`: `T(0,0)`
- `rank1`: `T(1,1)`, `T(1,0)`, `T(1,-1)`
- `rank2`: `T(2,2)`, `T(2,1)`, `T(2,0)`, `T(2,-1)`, `T(2,-2)`

## Inputs and checks

Pass all three arguments; use an empty numeric array such as `[]` for an omitted rank. Each input must be numeric. A nonempty `rank0`, `rank1`, or `rank2` must contain exactly `1`, `3`, or `5` elements, respectively. The source header describes these coefficients as vectors, while the implementation’s validation checks numeric type and element count, not vector shape.

## Cartesian matrix construction

The function initialises `M=zeros(3)` and adds each nonempty rank contribution. The following are the exact source coefficients and matrices:

```matlab
if ~isempty(rank0), M=M+rank0*eye(3); end
if ~isempty(rank1)
    M=M-(1/2)*[0 0 -1; 0 0 -1i; 1 1i 0]*rank1(1);
    M=M+(1/sqrt(2))*[0 -1i 0; 1i 0 0; 0 0 0]*rank1(2);
    M=M-(1/2)*[0 0 -1; 0 0 1i; 1 -1i 0]*rank1(3);
end
if ~isempty(rank2)
    M=M+(1/2)*[1 1i 0; 1i -1 0; 0 0 0]*rank2(1);
    M=M-(1/2)*[0 0 1; 0 0 1i; 1 1i 0]*rank2(2);
    M=M+(1/sqrt(6))*[-1 0 0; 0 -1 0; 0 0 2]*rank2(3);
    M=M+(1/2)*[0 0 1; 0 0 -1i; 1 -1i 0]*rank2(4);
    M=M+(1/2)*[1 -1i 0; -1i -1 0; 0 0 0]*rank2(5);
end
```

The output is a `3x3` matrix. The mapping uses the source’s factors `1/2`, `1/sqrt(2)`, and `1/sqrt(6)` and may produce complex matrix entries. No isotope-specific example appears in this source or its existing KB page.
