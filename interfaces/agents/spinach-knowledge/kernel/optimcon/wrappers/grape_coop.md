# kernel/optimcon/wrappers/grape_coop.m

- Signature: `[traj_data,fidelity,gradient]=grape_coop(phi_profile,spin_system)`

## Purpose

Pairs of cooperative pulses that may be used as components of a phase cycle. The pulses are designed to produce as much of the destination state as they can, and to have imputities of opposite sign. Adding the outcomes of the two experiments then destroys the impurities. Syntax: [traj_data,fidelity,gradient]=grape_coop(phi_profile,spin_system)

## Physical / mathematical content

The cooperative objective averages the requested primary transfer fidelities and subtracts the mean squared norm of the summed orthogonal impurities. Auxiliary impurity derivatives use real linear overlap, including when the primary fidelity is the absolute square of the overlap; target projection retains the target norm explicitly. Purely imaginary auxiliary overlaps and vanishing impurities have valid real-linear derivatives, including zero derivatives; the low-level engines return these without applying the optimiser's initial-gradient check.

- Optimal-control core routines. These files implement GRAPE-style objective evaluation, quasi-Newton search, line search, regularisation, distortion models, and waveform parameterisations.
- The control theory content is GRAPE: fidelity derivatives are propagated through a piecewise-constant pulse sequence so that waveform samples can be improved by gradient-based optimisation.

## Numerical / algorithmic content

## Parameters / inputs

- phi_profile -phase profiles of the two pulses,
- concatenated horizontally

## Outputs

- traj_data -trajectory information for both pulses
- fidelity -cooperative fidelity measure
- gradient -cooperative fidelity gradient
- Note: only phase-modulated point-to-point transformations are supported.
