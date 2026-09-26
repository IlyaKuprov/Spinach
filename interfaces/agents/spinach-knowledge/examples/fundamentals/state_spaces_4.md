# examples/fundamentals/state_spaces_4.m

- Signature: `state_spaces_4()`

## Purpose

Trajectory analysis for a MAS simulation of isotopically labelled glycine powder, starting from L+ on protons. Calculation time: hours, faster on a GPU.

## Physical / mathematical content

- Fundamentals examples. These are unit tests, convention checks, and pedagogical demonstrations of operator algebra, perturbation theory, tensor conventions, symmetry, quadrature, and numerical differentiation.
- Propagation is accelerated with a Krylov-subspace method, replacing direct matrix exponentiation by projection into a much smaller Arnoldi/Lanczos-type subspace.

## Numerical / algorithmic content

- A Krylov-subspace or Arnoldi construction is used to avoid forming or exponentiating very large dense propagators directly.

## Implementation structure

- Trajectory analysis for a MAS simulation of isotopically labelled
- glycine powder, starting from L+ on protons.
- Calculation time: hours, faster on a GPU.
- Spin system properties (PCM DFT calculation)
- Magnet field
- Basis set
- Force Krylov propagation
- This needs a GPU
- sys.enable={'gpu'};
- Spinach housekeeping
- Experiment setup
- Get the trajectory
