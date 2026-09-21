# examples/relaxation_theory/sle_solid_limit.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/relaxation_theory/sle_solid_limit.m`
- Signature: `sle_solid_limit()`
- Total lines: 73

## Purpose

Solid limit of Stochastic Liouville equation formalism. Calculation time: hours

## Physical / mathematical content

- Relaxation-theory examples. The mathematical backbone is Bloch-Redfield-Wangsness or stochastic Liouville theory, spectral densities, cross-correlation terms, motional models, and extraction of longitudinal/transverse decay behaviour from superoperators.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Implementation structure

- Solid limit of Stochastic Liouville equation formalism.
- Calculation time: hours
- Magnet field
- Isotopes
- Coupling Matrices
- Zeeman Interactions
- Basis set
- Spinach housekeeping
- SLE parameters
- Ranks and correlation times
- Start a figure
- Loop over correlation times

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `gauss2mhz()`, `create()`, `basis()`, `state()`, `kfigure()`, `scale_figure()`, `ranks()`, `tau_c()`, `gridfree()`, `subplot()`, `plot_1d()`, `ktitle()`, `num2str()`, `log10()`, `kxlabel()`.
