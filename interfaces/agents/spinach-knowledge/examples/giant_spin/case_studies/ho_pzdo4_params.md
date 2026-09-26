# examples/giant_spin/case_studies/ho_pzdo4_params.m

- Signature: `[ks,qs,bkq]=ho_pzdo4_params()`

## Purpose

Crystal field parameters of the Ho(pzdo)4 metal-organic framework in the extended Stevens operator convention, ranks 2 to 12, as computed at CASSCF level and used in the pulsed-field magnetisation case studies of (https://arxiv.org/abs/2609.16352).

```matlab
[ks,qs,bkq]=ho_pzdo4_params()
```

## Outputs

- ks -row of Stevens operator ranks
- qs -row of Stevens operator projections
- bkq -row of coefficients, cm^-1

## Physical / mathematical content

- 90 coefficients: even ranks 2, 4, 6, 8, 10, and 12 with every projection q from -k to k, in cm^-1, taken from the input file of the paper. The odd-q and high-rank values are numerically tiny and are kept as computed.
- Both Ho scripts assemble, for each rank k, a vector of length 2k+1 indexed by q+k+1 from these rows, convert it with `icm2hz` and `stev2sph`, and place the result in `inter.giant.coeff{1}{k}` with zero Euler angles.
