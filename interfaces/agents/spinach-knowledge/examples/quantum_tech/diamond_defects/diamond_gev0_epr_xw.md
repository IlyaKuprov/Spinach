# examples/quantum_tech/diamond_defects/diamond_gev0_epr_xw.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/quantum_tech/diamond_defects/diamond_gev0_epr_xw.m`
- Signature: `diamond_gev0_epr_xw()`
- Total lines: 67

## Purpose

Field-swept powder EPR spectra of GeV0 centre in diamond at X and W bands. Calculation time: seconds.

## Physical / mathematical content

- Quantum-technology examples. The files in this area model cavity QED, transmon qubits, NV centres, and related effective Hamiltonians. The recurring mathematics is finite-dimensional quantum dynamics with ladder operators, rotating-wave-style couplings, anharmonic oscillator terms, avoided crossings, and coherent control in coupled few-mode systems.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Implementation structure

- Field-swept powder EPR spectra of GeV0 centre
- in diamond at X and W bands.
- Calculation time: seconds.
- Set GeV0 centre model parameters.
- Build the spin system.
- Field sweep
- Define the basis set
- Run Spinach housekeeping
- Set common EPR parameters
- Set X-band parameters
- Run the X-band simulation
- Plot the X-band spectrum

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `diamond_gev0()`, `create()`, `basis()`, `fieldsweep()`, `kfigure()`, `scale_figure()`, `subplot()`, `kxlabel()`, `kylabel()`, `ktitle()`.
