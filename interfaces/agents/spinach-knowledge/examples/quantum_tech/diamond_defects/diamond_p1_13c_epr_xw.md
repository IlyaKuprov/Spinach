# examples/quantum_tech/diamond_defects/diamond_p1_13c_epr_xw.m

- Signature: `diamond_p1_13c_epr_xw()`

## Purpose

Field-swept powder EPR spectra of a P1 centre in 13C-enriched diamond at X and W bands. Calculation time: minutes.

## Physical / mathematical content

- Quantum-technology examples. The files in this area model cavity QED, transmon qubits, NV centres, and related effective Hamiltonians. The recurring mathematics is finite-dimensional quantum dynamics with ladder operators, rotating-wave-style couplings, anharmonic oscillator terms, avoided crossings, and coherent control in coupled few-mode systems.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Implementation structure

- Field-swept powder EPR spectra of a P1 centre
- in 13C-enriched diamond at X and W bands.
- Calculation time: minutes.
- Set P1 model parameters.
- Build the spin system.
- Field sweep
- Define the basis set
- Run Spinach housekeeping
- Leave only 14N nucleus and 13C nuclei with hyperfines larger than
- a_iso > 8 MHz
- Set common EPR parameters
- Set X-band parameters
