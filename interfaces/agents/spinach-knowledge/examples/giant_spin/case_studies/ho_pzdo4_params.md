# examples/giant_spin/case_studies/ho_pzdo4_params.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/giant_spin/case_studies/ho_pzdo4_params.m`
- Signature: `[ks,qs,bkq]=ho_pzdo4_params()`
- Total lines: 36

## Purpose

Crystal field parameters of the Ho(pzdo)4 metal-organic framework in the extended Stevens operator convention, ranks 2 to 12, as computed at CASSCF level and used in the pulsed-field magnetisation case studies of https://arxiv.org/abs/2609.16352 (their input file).

## Physical / mathematical content

- Giant-spin examples. The effective model treats lanthanides or high-spin centres using crystal-field / Stevens-operator Hamiltonians, Zeeman splitting, and magnetisation dynamics.

## Numerical / algorithmic content


## Code-derived implementation details

### Comment-guided execution stages

- Lines 21-24: Crystal field parameters, cm^-1, ranks 2 to 12 in Stevens operator convention; implemented by `ks=[2 2 2 2 2 4 4 4 4 4 4 4 4 4 6 6 6 6 6 6 6 6 6 6 6 6 6 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 12 12 12 12 12…`.

### Key state/data transformations

- Lines 22-24: computes `ks` using `ks=[2 2 2 2 2 4 4 4 4 4 4 4 4 4 6 6 6 6 6 6 6 6 6 6 6 6 6 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 12 12 12 12 12…`.
- Lines 25-27: computes `qs` using `qs=[-2 -1 0 1 2 -4 -3 -2 -1 0 1 2 3 4 -6 -5 -4 -3 -2 -1 0 1 2 3 4 5 6 -8 -7 -6 -5 -4 -3 -2 -1 0 1 2 3 4 5 6 7 8 -10 -9 -8 -7 -6 -5 -4 -3 -2 -1 0 1 2 3 4 5 6 7 8 9 10 -12…`.
- Lines 28-33: computes `bkq` using `bkq=[3.4994750592E-02 -2.7086096290E-06 2.7144135449E-02 -3.8607263323E-06 8.5970155653E-02 4.9292781766E-03 3.7420773069E-07 -6.3647621451E-04 -5.8482618912E-07 2.58411…`.

## Syntax

```matlab
[ks,qs,bkq]=ho_pzdo4_params()
```

## Outputs

- ks -row of Stevens operator ranks
- qs -row of Stevens operator projections
- bkq -row of coefficients, cm^-1

## Implementation structure

- Crystal field parameters of the Ho(pzdo)4 metal-organic framework
- in the extended Stevens operator convention, ranks 2 to 12, as
- computed at CASSCF level and used in the pulsed-field magnetisation
- case studies of https://arxiv.org/abs/2609.16352 (their input file).
- [ks,qs,bkq]=ho_pzdo4_params()
- ks -row of Stevens operator ranks
- qs -row of Stevens operator projections
- bkq -row of coefficients, cm^-1
- Crystal field parameters, cm^-1, ranks 2 to 12 in Stevens operator convention
