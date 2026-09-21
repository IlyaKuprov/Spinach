# examples/singlet_states/singlet_imaging_1.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/singlet_states/singlet_imaging_1.m`
- Signature: `singlet_imaging_1()`
- Total lines: 160

## Purpose

Singlet imaging in a system with one-dimensional diffusion and flow. Calculation time: minutes

## Physical / mathematical content

- Long-lived singlet-state examples. The central concept is symmetry-protected or nearly symmetry-protected two-spin order that relaxes much more slowly than ordinary Zeeman magnetisation. Files here often analyse singlet-triplet subspaces, state conversion sequences, and relaxation leakage channels.
- The relaxation model is Redfield-type perturbation theory: fluctuating interactions enter through correlation functions or spectral densities and generate a linear relaxation superoperator.
- The relevant state manifold is the singlet/triplet decomposition, where permutation symmetry controls selection rules, relaxation susceptibility, and convertibility to ordinary magnetisation.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- The file also defines local helper function(s): `tube_flow()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Implementation structure

- Singlet imaging in a system with one-dimensional
- diffusion and flow.
- Calculation time: minutes
- Spin system and interactions
- Relaxation theory
- Relaxation superoperator accuracy
- Algorithmic options
- Basis set
- Spinach housekeeping
- Sample geometry
- Sequence parameters
- Assumptions

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `assume()`, `relaxation()`, `tube()`, `state()`, `tube_flow()`, `operator()`, `speye()`, `pulse_shape()`, `shaped_pulse_af()`, `step()`, `evolution()`, `imaging()`, `singlet()`, `kfigure()`.
