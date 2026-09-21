# examples/relaxation_theory/trosy_methyl.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/relaxation_theory/trosy_methyl.m`
- Signature: `trosy_methyl()`
- Total lines: 155

## Purpose

Methyl trosy in a rapidly rotating 13CH3 group of a slowly tumbling protein, simulated using the Fokker-Planck formalism.

## Physical / mathematical content

- Relaxation-theory examples. The mathematical backbone is Bloch-Redfield-Wangsness or stochastic Liouville theory, spectral densities, cross-correlation terms, motional models, and extraction of longitudinal/transverse decay behaviour from superoperators.
- The file uses a Fokker-Planck-style enlarged state space in which spatial or orientational coordinates are promoted to extra dimensions and coupled to spin dynamics through differential operators.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Implementation structure

- Methyl trosy in a rapidly rotating 13CH3 group
- of a slowly tumbling protein, simulated using
- the Fokker-Planck formalism.
- Magnet field
- Cartesian coordinates
- Absolute shielding tensors (DFT)
- Convert shielding tensors into
- chemical shift tensors and put
- them on resonance
- Methyl proton chemical shifts (guess)
- J-couplings
- Spin system instances for the three methyl rotamers

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `remtrace()`, `me_xyz()`, `create()`, `basis()`, `kfigure()`, `scale_figure()`, `state()`, `gridfree()`, `subplot()`, `plot_1d()`.
