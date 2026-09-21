# examples/singlet_states/eigenstate_analysis.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/singlet_states/eigenstate_analysis.m`
- Signature: `eigenstate_analysis()`
- Total lines: 100

## Purpose

Stationary state analysis for the spin system of allyl pyruvate, finding out which component of the singlet state commutes with the drift Hamiltonian.

## Physical / mathematical content

- Long-lived singlet-state examples. The central concept is symmetry-protected or nearly symmetry-protected two-spin order that relaxes much more slowly than ordinary Zeeman magnetisation. Files here often analyse singlet-triplet subspaces, state conversion sequences, and relaxation leakage channels.
- The relevant state manifold is the singlet/triplet decomposition, where permutation symmetry controls selection rules, relaxation susceptibility, and convertibility to ordinary magnetisation.

## Numerical / algorithmic content

- An eigenvalue problem is solved or analysed, so the file is extracting spectra, stationary states, avoided crossings, or modal structure from the effective Hamiltonian or superoperator.

## Implementation structure

- Stationary state analysis for the spin system of allyl pyruvate,
- finding out which component of the singlet state commutes with
- the drift Hamiltonian.
- Get the spin system from Anu's fits
- Set the magnet
- Spinach housekeeping
- Pick out the required 13C isotopomer
- Generate the basis
- Get isotropic Hamiltonian
- Tidy up rounding errors
- Get the singlet state
- Report the norm

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `allyl_pyruvate()`, `create()`, `dilute()`, `basis()`, `assume()`, `hamiltonian()`, `ctranspose()`, `singlet()`, `report()`, `num2str()`, `remncomm()`, `remtrace()`, `state()`, `frqoffset()`, `operator()`.
