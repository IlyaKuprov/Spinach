# examples/singlet_states/decoherence_diacetylene.m

- Signature: `decoherence_diacetylene()`

## Purpose

Long-lived spin states in the diacetylene molecule (2 protons, 4 carbons, 4096-dimensional Liouville space). The relaxation superoperator accounts for every dipolar coupling and every CSA tensor in the system. Calculation time: seconds

## Physical / mathematical content

For diacetylene with 2 protons and 4 carbons, the Redfield relaxation superoperator includes every dipolar coupling and CSA tensor, and the calculation examines the singlet between the two central carbons.

## Numerical / algorithmic content

In 4,096-dimensional Liouville space, the code lists 20 smallest-magnitude relaxation eigenvalues, evaluates the normalized singlet self-relaxation rate, and analyzes two slowly relaxing eigenvectors.

## Implementation structure

The code imports vacuum-DFT spin data, sets the field to 14.1 T, uses Redfield relaxation with zero equilibrium, lab-frame retention and tau_c=100e-12, sets the relaxation integration and zero tolerances to 1e-5, and builds an untruncated sphten-liouv basis.
