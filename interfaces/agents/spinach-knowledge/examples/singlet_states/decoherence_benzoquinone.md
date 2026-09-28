# examples/singlet_states/decoherence_benzoquinone.m

- Signature: `decoherence_benzoquinone()`

## Purpose

Long-lived spin states in the para-benzoquinone molecule (4 protons, 256-dimensional Liouville space). The relaxation superoperator accounts for every dipolar coupling and every CSA tensor in the system. Calculation time: seconds

## Physical / mathematical content

The example examines long-lived spin states in para-benzoquinone’s four-proton system using a Redfield relaxation superoperator that includes every dipolar coupling and CSA tensor.

## Numerical / algorithmic content

It uses a 256-dimensional Liouville space, a 1.0 T field, a 100 ps correlation time, 1e-5 relaxation tolerances, and lists the 20 smallest-magnitude relaxation eigenvalues.

## Implementation structure

The function imports vacuum-DFT spin data, constructs a complete spherical-tensor Liouville basis, builds the relaxation superoperator, and extracts its eigenvalues.
