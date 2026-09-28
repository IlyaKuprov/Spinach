# examples/relaxation_theory/nz_vs_redfield_1.m

- Signature: `nz_vs_redfield_1()`

## Purpose

Compares Redfield relaxation with on-shell and off-shell Nakajima–Zwanzig (NZ) kernels for a two-spin system with dipolar and CSA cross-correlations. It evaluates the on-shell/Redfield difference, the off-shell/Redfield difference on the zero-frequency subspace and across the full superoperator, the dependence on correlation time, and the effect of an NZ lifetime shift. The source estimates minutes of calculation time.

## Physical / mathematical content

The system is 1H/13C at 14.1 T, with shielding principal values `[7 15 -22]` and `[11 18 -29]`, Euler angles as specified in the script, and a 1.02-coordinate-unit separation for the dipolar coupling. Common settings are zero equilibrium, lab-frame retention, and `rlx_dfs='keep'`; the complete `sphten-liouv` basis is used without approximation.

## Numerical / algorithmic content

Redfield and NZ relaxation superoperators are converted to full matrices. At zero shift the script compares on-shell NZ with Redfield; for off-shell NZ it also compares their action on the nullspace of the lab-frame coherent Hamiltonian. It scans correlation times `[2.5, 5, 10, 20]` ps and lifetime shifts `[0, 1e9, 1e10, 1e11]` Hz, recording relative matrix-norm differences and the largest absolute diagonal relaxation rate. These computed trends are plotted; no numerical outcomes are hard-coded in the page.

## Implementation structure

The code builds separate spin systems and relaxation matrices for Redfield, on-shell NZ (`nz_onshell=true`, `nz_shift=0`), and off-shell NZ. It computes the relative norm comparisons, scans the correlation-time and lifetime-shift grids, and plots the two trends in side-by-side panels.
