# examples/nmr_solids/dor_powder_nav_fplanck_freq.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_solids/dor_powder_nav_fplanck_freq.m`
- Signature: `dor_powder_nav_fplanck_freq()`
- Total lines: 63

## Purpose

Double angle spinning spectrum of N-acetylvaline 14N nucleus using 1D Fokker-Planck equation and a spherical grid. The cal- culation includes the second-order quadrupolar shift and the third-order lineshape. Frequency-domain detection within the user-specified frequency interval. Note: slower spinning rates and larger NQIs require larger ranks and spherical grids. At the moment the spinning frequencies are set artifi

## Physical / mathematical content

- Solid-state NMR examples. The key physics is anisotropic spin interactions under static or magic-angle-spinning conditions: chemical-shift anisotropy, dipolar coupling, quadrupolar coupling, cross-polarisation, and orientation averaging using Floquet, Fokker-Planck, or direct powder quadrature formalisms.
- The file uses a Fokker-Planck-style enlarged state space in which spatial or orientational coordinates are promoted to extra dimensions and coupled to spin dynamics through differential operators.
- Quadrupolar physics is relevant: nuclei with spin > 1/2 interact with the electric field gradient tensor, introducing second-rank anisotropy, asymmetry, and overtone or MQ phenomena.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Implementation structure

- Double angle spinning spectrum of N-acetylvaline 14N nucleus
- using 1D Fokker-Planck equation and a spherical grid. The cal-
- culation includes the second-order quadrupolar shift and the
- third-order lineshape. Frequency-domain detection within the
- user-specified frequency interval.
- Note: slower spinning rates and larger NQIs require larger
- ranks and spherical grids. At the moment the spinning
- frequencies are set artificially too high to reduce
- the simulatio ntime in this example.
- Calculation time: minutes
- System specification
- Relaxation theory

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `eeqq2nqi()`, `create()`, `basis()`, `state()`, `doublerot()`, `kfigure()`, `plot_1d()`.
