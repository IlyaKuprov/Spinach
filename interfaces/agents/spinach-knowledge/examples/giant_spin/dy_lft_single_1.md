# examples/giant_spin/dy_lft_single_1.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/giant_spin/dy_lft_single_1.m`
- Signature: `dy_lft_single_1()`
- Total lines: 97

## Purpose

Reproduction of MOLCAS results with the Ligand Field Theory model for a single Dy(III) ion. Calculation time: seconds

## Physical / mathematical content

- Giant-spin examples. The effective model treats lanthanides or high-spin centres using crystal-field / Stevens-operator Hamiltonians, Zeeman splitting, and magnetisation dynamics.

## Numerical / algorithmic content

- An eigenvalue problem is solved or analysed, so the file is extracting spectra, stationary states, avoided crossings, or modal structure from the effective Hamiltonian or superoperator.

## Implementation structure

- Reproduction of MOLCAS results with the Ligand Field Theory model
- for a single Dy(III) ion.
- Calculation time: seconds
- Magnetic field
- Single Dy ion
- Real g-tensor
- Rotate the ligand field into the molecular frame
- Liza -this needs more decimal places
- Ligand field parameters (MOLCAS)
- Convert to irreducible spherical tensors
- Supply to Spinach
- Formalism specification

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `dcm2euler()`, `icm2hz()`, `wigner()`, `stev2sph()`, `create()`, `basis()`, `geffect()`, `num2str()`.
