# examples/optimal_control/case_studies/Tosner_JMR_2009/coherence_transfer.m

- Signature: `coherence_transfer()`

## Purpose

This first optimal-control example, associated with http://dx.doi.org/10.1016/j.jmr.2008.11.020, uses an on-resonance heteronuclear two-spin system (1H–13C) with a 140 Hz scalar coupling. The goal is to transfer transverse magnetisation from proton to carbon (Hx → Cx) over a fixed period T = 1/J.

## Physical / mathematical content

- The initial and target states are normalised proton and carbon Lx states, respectively. The controls use Lx and Ly operators on each nucleus, alongside the NMR drift Hamiltonian.
- The pulse spans 150 equal time steps totalling 1/J. The control structure specifies `lbfgs` as the optimiser, a maximum of 200 iterations, power levels `2*pi*linspace(10,1000,10)`, and an `NS` penalty with weight 0.01. Optimisation calls `fmaxnewton` with `@grape_xy`.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Implementation structure

- Set the magnetic field to 14.1 T, both chemical shifts to 0 ppm, and the 1H–13C scalar coupling to 140 Hz.
- Create the spin system with the `sphten-liouv` formalism and `none` approximation; construct and normalise the proton Lx initial state and carbon Lx target state.
- Configure the drift and four control operators, generate a random 4-by-150 initial guess divided by 10, and run the optimisation. Requested plots are `xy_controls`, `spectrogram`, and `robustness`.
