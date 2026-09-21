# examples/relaxation_theory/cpmg_echo_train.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/relaxation_theory/cpmg_echo_train.m`
- Signature: `cpmg_echo_train()`
- Total lines: 55

## Purpose

CPMG echo train in a powder. Calculation time: seconds

## Physical / mathematical content

- Relaxation-theory examples. The mathematical backbone is Bloch-Redfield-Wangsness or stochastic Liouville theory, spectral densities, cross-correlation terms, motional models, and extraction of longitudinal/transverse decay behaviour from superoperators.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Implementation structure

- CPMG echo train in a powder.
- Calculation time: seconds
- System specification
- Basis set
- Relaxation theory
- Algorithmic options
- Spinach housekeeping
- Experiment setup
- Simulation
- Plotting

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `state()`, `operator()`, `powder()`, `kfigure()`, `scale_figure()`, `kxlabel()`, `kylabel()`.
