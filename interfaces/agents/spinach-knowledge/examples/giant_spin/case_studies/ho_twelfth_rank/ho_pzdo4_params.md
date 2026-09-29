# examples/giant_spin/case_studies/ho_twelfth_rank/ho_pzdo4_params.m

- MATLAB implementation: [examples/giant_spin/case_studies/ho_twelfth_rank/ho_pzdo4_params.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/giant_spin/case_studies/ho_twelfth_rank/ho_pzdo4_params.m)

- Signature: `[ks,qs,bkq]=ho_pzdo4_params()`

## Purpose

Returns the CASSCF crystal-field coefficients for the Ho(pzdo)4 metal-organic framework, used in the pulsed-field magnetisation examples associated with [arXiv:2609.16352](https://arxiv.org/abs/2609.16352). The coefficients use the extended Stevens operator convention and include even ranks 2, 4, 6, 8, 10, and 12.

## Interface and physical content

`ks`, `qs`, and `bkq` are row vectors: rank, projection, and coefficient, respectively. For each rank k, the source lists every q from -k through k, for 90 coefficients in total. The coefficients `bkq` are in cm^-1. Numerically small odd-q and high-rank values are retained as supplied; the function does not threshold or recalculate them.

In the two Ho simulations, the rows are regrouped rank by rank into vectors of length 2k+1, with element q+k+1 corresponding to projection q. Each rank's vector is converted first with `icm2hz`, then with `stev2sph`, and stored as `inter.giant.coeff{1}{k}`; the associated Euler angles are [0 0 0].

## Output and limits

The function returns the parameter rows only. It does not build a spin system, propagate magnetisation, or produce a plotted or numerical simulation result. The cited paper is the source context for the CASSCF parameter set: https://arxiv.org/abs/2609.16352.
