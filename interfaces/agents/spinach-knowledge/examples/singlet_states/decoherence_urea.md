# examples/singlet_states/decoherence_urea.m

- Signature: `decoherence_urea()`

## Purpose

A demonstration that the nitrogen singlet state in urea is not long-lived. The relaxation superoperator accounts for every di- polar coupling and every CSA tensor in the system. Calculation time: seconds

## Physical / mathematical content

The example probes relaxation of the 15N singlet state in urea, including all dipolar couplings and CSA tensors in the system.

## Numerical / algorithmic content

At 1.0 T, it builds a Redfield relaxation superoperator with a 100 ps correlation time, zero equilibrium and lab-frame terms; the relaxation integration and zero tolerances are 1e-5.

## Implementation structure

It reads urea coordinates, shifts, J-couplings and CSAs from `../standard_systems/urea.log`, uses a complete sphten-liouv basis, and prints `norm(R*Lz)/norm(Lz)` and `norm(R*S)/norm(S)` for 15N Lz and the singlet of spins 1 and 4.
