# examples/singlet_states/singlet_imaging_1.m

- Signature: `singlet_imaging_1()`

## Purpose

Singlet imaging in a system with one-dimensional diffusion and flow. Calculation time: minutes

## Physical / mathematical content

- The model uses two 13C spins with scalar coupling 55 Hz and opposite Zeeman offsets, 0.03 and -0.03, under a Redfield relaxation model with zero equilibrium, secular retention, and correlation time `1e-9`. The example compares singlet-state and transverse-magnetisation imaging trajectories.
- The 2D sample grid is `150 x 15` over dimensions `[0.10 0.015]`; flow is along the first coordinate at `-6e-2`, and diffusion is `3.6e-6` along that coordinate only.

## Numerical / algorithmic content

- The script builds a spatial Liouvillian from the Hamiltonian, flow, relaxation, and kinetic terms, then applies a shaped pulse and an M2S sequence to prepare singlet order. `imaging` propagates both singlet and magnetisation states for display as separate tube images.

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
