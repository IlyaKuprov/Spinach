# examples/singlet_states/decoherence_naphthalenetetrone.m

- Signature: `decoherence_naphthalenetetrone()`

## Purpose

Long-lived spin states in the naphthalenetetrone molecule. (4 protons, 256-dimensional Liouville space). The relaxation superoperator accounts for every dipolar coupling and every CSA tensor in the system. Calculation time: seconds

## Physical / mathematical content

The four-proton naphthalenetetrone system includes every dipolar coupling and CSA tensor in a Redfield relaxation model at 1.0 T, with a 100 ps correlation time.

## Numerical / algorithmic content

In the complete 256-dimensional Liouville space, `eigs` computes the twenty smallest-magnitude eigenvalues of `R-speye(size(R))`, then adds 1 to report relaxation rates in Hz.

## Implementation structure

The function imports vacuum-DFT coordinates, shifts, J-couplings and CSAs from `../standard_systems/naphthalenetetrone.log`, selects 1H spins, sets zero equilibrium and lab-frame relaxation with both relaxation tolerances at 1e-5, and builds the sphten-liouv basis without approximation.
