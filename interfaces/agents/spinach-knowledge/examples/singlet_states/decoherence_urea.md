# examples/singlet_states/decoherence_urea.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/singlet_states/decoherence_urea.m`
- Signature: `decoherence_urea()`
- Total lines: 52

## Purpose

A demonstration that the nitrogen singlet state in urea is not long-lived. The relaxation superoperator accounts for every di- polar coupling and every CSA tensor in the system. Calculation time: seconds

## Physical / mathematical content

- Long-lived singlet-state examples. The central concept is symmetry-protected or nearly symmetry-protected two-spin order that relaxes much more slowly than ordinary Zeeman magnetisation. Files here often analyse singlet-triplet subspaces, state conversion sequences, and relaxation leakage channels.
- The relaxation model is Redfield-type perturbation theory: fluctuating interactions enter through correlation functions or spectral densities and generate a linear relaxation superoperator.
- Chemical-shift anisotropy is present: shielding is treated as a second-rank tensor whose orientation relative to the field or rotor axis modulates line shapes and transfer dynamics.
- The relevant state manifold is the singlet/triplet decomposition, where permutation symmetry controls selection rules, relaxation susceptibility, and convertibility to ordinary magnetisation.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Implementation structure

- A demonstration that the nitrogen singlet state in urea is not
- long-lived. The relaxation superoperator accounts for every di-
- polar coupling and every CSA tensor in the system.
- Calculation time: seconds
- Read the spin system (coordinates, chemical shifts,
- J-couplings and CSAs) from a vacuum DFT calculation
- Set magnet field to 1.0 Tesla
- Tighten up the tolerances
- Set relaxation theory parameters
- Relaxation superoperator accuracy
- Use complete basis set
- Spinach housekeeping

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `g2spinach()`, `gparse()`, `create()`, `basis()`, `relaxation()`, `state()`, `singlet()`, `num2str()`.
