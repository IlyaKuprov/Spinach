# examples/relaxation_theory/nz_vs_redfield_1.m

[Source file](https://github.com/IlyaKuprov/Spinach/blob/main/examples/relaxation_theory/nz_vs_redfield_1.m) · Signature: `nz_vs_redfield_1()`

## Purpose and model

This example compares Redfield relaxation with two Nakajima–Zwanzig (NZ) kernels for a two-spin system with dipolar and chemical-shift-anisotropy cross-correlations. The source distinguishes the off-shell resolvent kernel from the on-shell back-rotated kernel. Its header states the expected relationships: at zero shift the on-shell kernel reproduces Redfield theory; the off-shell kernel agrees with Redfield on the zero-frequency subspace and differs to first order in omega times correlation time on coherences; a lifetime shift suppresses rates by moving the kernel off the real axis. These are source-described theoretical comparisons, not experimental observations.

The model uses `1H` and `13C` with magnet setting 14.1. The shielding principal-value inputs are [7, 15, -22] and [11, 18, -29], with Euler-angle inputs [pi/3, pi/4, pi/5] and [pi/6, pi/7, pi/8]. Coordinates are (0,0,0) and (0,0,1.02); the source does not annotate a unit for these entries. Common relaxation settings are zero equilibrium, lab-frame retention, and keeping the density-fluctuation settings. The Redfield and NZ base correlation-time inputs are 200e-12. The source explicitly labels seconds on the correlation-time plot axis; it does not otherwise attach a unit to the base input in the parameter assignment.

## Comparisons and displayed output

The script constructs full relaxation superoperators for Redfield, on-shell NZ with `nz_onshell=true` and zero shift, and off-shell NZ. It prints relative 1-norm differences for on-shell versus Redfield, off-shell versus Redfield on the zero-frequency (kite) subspace, and off-shell versus Redfield over the full superoperator. It then scans correlation times 2.5e-12, 5e-12, 10e-12, and 20e-12, comparing off-shell NZ and Redfield at each value. A second scan applies lifetime shifts 0, 1e9, 1e10, and 1e11 Hz and records the largest absolute diagonal relaxation rate. Two panels plot the calculated relative difference versus correlation time and largest rate versus lifetime shift. The source gives no measured data and the example page reports no run-specific numerical output; the plotted curves and printed values are calculations produced when the script runs.
