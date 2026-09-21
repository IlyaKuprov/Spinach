# examples/giant_spin/nuclear_relaxation_1.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/giant_spin/nuclear_relaxation_1.m`
- Signature: `nuclear_relaxation_1()`
- Total lines: 148

## Purpose

Nuclear relaxation rates using the adiabatic elimination method for a rapidly relaxing Dy(III) ion with a user-specified ZFS. Calculation time: minutes

## Physical / mathematical content

- Giant-spin examples. The effective model treats lanthanides or high-spin centres using crystal-field / Stevens-operator Hamiltonians, Zeeman splitting, and magnetisation dynamics.

## Numerical / algorithmic content

- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.

## Implementation structure

- Nuclear relaxation rates using the adiabatic elimination method
- for a rapidly relaxing Dy(III) ion with a user-specified ZFS.
- Calculation time: minutes
- Magnetic field
- Dy(III) ion and a proton
- Electron g-tensor
- Spin-orbit corrections
- to the DD couplings
- Nuclear shift tensor
- Rotate the ligand field into the molecular frame
- Liza -this needs more decimal places
- Ligand field parameters (MOLCAS)

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `dcm2euler()`, `icm2hz()`, `wigner()`, `stev2sph()`, `create()`, `basis()`, `hamiltonian()`, `assume()`, `relaxation()`, `load()`, `orientation()`, `alphas()`, `betas()`, `gammas()`, `adelim()`, `weights()`.
