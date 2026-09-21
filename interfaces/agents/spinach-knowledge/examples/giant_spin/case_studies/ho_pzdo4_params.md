# examples/giant_spin/case_studies/ho_pzdo4_params.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/giant_spin/case_studies/ho_pzdo4_params.m`
- Signature: `[ks,qs,bkq]=ho_pzdo4_params()`
- Total lines: 53

## Purpose

Crystal field parameters of the Ho(pzdo)4 metal-organic framework in the extended Stevens operator convention, ranks 2 to 12, as computed at CASSCF level and used in the pulsed-field magnetisation case studies of (https://arxiv.org/abs/2609.16352).

## Physical / mathematical content

- Giant-spin examples. The effective model treats lanthanides or high-spin centres using crystal-field / Stevens-operator Hamiltonians, Zeeman splitting, and magnetisation dynamics.

## Numerical / algorithmic content


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
- case studies of (https://arxiv.org/abs/2609.16352).
- [ks,qs,bkq]=ho_pzdo4_params()
- ks -row of Stevens operator ranks
- qs -row of Stevens operator projections
- bkq -row of coefficients, cm^-1
- Crystal field parameters, cm^-1, ranks 2 to 12 in Stevens operator convention
