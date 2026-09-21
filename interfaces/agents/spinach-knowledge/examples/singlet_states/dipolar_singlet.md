# examples/singlet_states/dipolar_singlet.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/singlet_states/dipolar_singlet.m`
- Signature: `dipolar_singlet()`
- Total lines: 48

## Purpose

A demonstration that the two-spin singet state is immune to dipolar relaxation. Full Redfield superoperator for dipolar relaxation in liquid state is computed and the norm of its action on a singlet state is printed to the console. Calculation time: seconds

## Physical / mathematical content

- Long-lived singlet-state examples. The central concept is symmetry-protected or nearly symmetry-protected two-spin order that relaxes much more slowly than ordinary Zeeman magnetisation. Files here often analyse singlet-triplet subspaces, state conversion sequences, and relaxation leakage channels.
- The relaxation model is Redfield-type perturbation theory: fluctuating interactions enter through correlation functions or spectral densities and generate a linear relaxation superoperator.
- The spin physics includes through-space magnetic dipole-dipole coupling, a rank-2 anisotropic interaction with strong orientation dependence and characteristic secular/non-secular structure.
- The relevant state manifold is the singlet/triplet decomposition, where permutation symmetry controls selection rules, relaxation susceptibility, and convertibility to ordinary magnetisation.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Implementation structure

- A demonstration that the two-spin singet state is immune
- to dipolar relaxation. Full Redfield superoperator for
- dipolar relaxation in liquid state is computed and the
- norm of its action on a singlet state is printed to the
- console.
- Calculation time: seconds
- System specification
- Relaxation theory parameters
- Relaxation superoperator accuracy
- Proximity cut-off
- Basis set
- Spinach housekeeping

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `relaxation()`, `singlet()`, `num2str()`.
