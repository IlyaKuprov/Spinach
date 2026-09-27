# examples/quantum_tech/optomechanics/optomech_sideband.m

- Signature: `optomech_sideband()`

## Purpose

Optomechanical sideband transfer from a phonon Fock state to a driven cavity. A red-detuned coherent cavity drive activates the beam-splitter component of radiation-pressure coupling, exchanging the mechanical mode's second Fock state with the cavity field. The model and parameters come from the propagation test set of QuantumPropagators.jl; its dimensionless quantities are divided by 2π to cancel the Hz convention. Calculation time: seconds.

## Physical / mathematical content

- The model has a five-level cavity and an eleven-level phonon mode. Radiation-pressure coupling links cavity photon number to the phonon coordinate; the coherent cavity drive produces the sideband exchange.

## Numerical / algorithmic content

- In the `zeeman-hilb` formalism with no basis approximation, the example builds the Hamiltonian, forms a one-step propagator, and iterates it over 250 steps of 0.2 dimensionless time units. It tracks both mode occupations and checks trace preservation, the initial phonon occupation, and transfer into the cavity.

## Implementation structure

- Both mode frequencies are `10/(2*pi)`; the longitudinal coupling is `-sqrt(2)/(2*pi)`. The drive adds `2*(C+A)` for cavity mode 1. The initial state is `{'BL1','BL3'}` on modes 1 and 2 (empty cavity and second phonon Fock state). The plotted observables are the cavity and mechanical-mode occupations.
