# examples/giant_spin/triple_dy_magn.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/giant_spin/triple_dy_magn.m`
- Signature: `triple_dy_magn()`
- Total lines: 133

## Purpose

Simulation of a finite-speed magnetic field sweep experiment for a single crystal of a triple-Dy triangular complex in a micro-SQUID, see Figure S24 in the Supplementary Information of Ligand field parameters and g-tensor for the J=15/2 ground term were computed using the SINGLE_ANISO routine in MOLCAS. Calculation time: hours

## Physical / mathematical content

- Giant-spin examples. The effective model treats lanthanides or high-spin centres using crystal-field / Stevens-operator Hamiltonians, Zeeman splitting, and magnetisation dynamics.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Implementation structure

- Simulation of a finite-speed magnetic field sweep experiment for a
- single crystal of a triple-Dy triangular complex in a micro-SQUID,
- see Figure S24 in the Supplementary Information of
- Ligand field parameters and g-tensor for the J=15/2 ground term were
- computed using the SINGLE_ANISO routine in MOLCAS.
- Calculation time: hours
- Three J=15/2 dysprosium atoms
- g-tensor eigenvalues
- Spin-orbit corrections
- to the DD couplings
- g-tensor eigenvectors
- g-tensor matrix

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `euler2dcm()`, `icm2hz()`, `dcm2euler()`, `wigner()`, `stev2sph()`, `create()`, `basis()`, `fieldscan_magn()`, `kfigure()`, `kxlabel()`, `kylabel()`.
