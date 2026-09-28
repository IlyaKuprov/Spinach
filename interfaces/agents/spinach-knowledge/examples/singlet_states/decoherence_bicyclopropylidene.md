# examples/singlet_states/decoherence_bicyclopropylidene.m

- Signature: `decoherence_bicyclopropylidene()`

## Purpose

Long-lived spin states in the bicyclopropylidene molecule (8 protons, 65536-dimensional Liouville space). The relaxation superoperator accounts for every dipolar coupling and every CSA tensor in the system. Calculation time: hours

## Physical / mathematical content

The eight-proton bicyclopropylidene model includes every dipolar coupling and CSA tensor in its relaxation superoperator.

## Numerical / algorithmic content

Redfield relaxation with a 100 ps correlation time is evaluated in a complete 65,536-dimensional Liouville space, and the twenty smallest-magnitude eigenvalues are listed as relaxation rates in Hz.

## Implementation structure

The function imports vacuum-DFT coordinates, shifts, J-couplings and CSAs, sets a 1.0 T field, zero equilibrium and lab-frame relaxation, and uses 1e-5 integration and zero tolerances before building the relaxation superoperator.
