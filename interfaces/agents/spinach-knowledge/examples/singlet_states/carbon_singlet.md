# examples/singlet_states/carbon_singlet.m

- Signature: `carbon_singlet()`

## Purpose

Calculate the singlet relaxation rate for the two triple-bond carbons in cis-dimethylbut-2-ynedioate using magnetic parameters computed with DFT. Calculation time: seconds.

## Physical / mathematical content

- The system contains two `13C` spins in a 14.1 T magnetic field, with specified Zeeman matrices and coordinates.
- Relaxation uses the Redfield model with a 100 ps correlation time, zero equilibrium, and lab-frame terms retained.
- The calculation evaluates relaxation-superoperator matrix elements for normalized longitudinal magnetization, `<Sz|R|Sz>`, and the normalized singlet state, `<singlet|R|singlet>`.

## Numerical / algorithmic content

- The basis uses `sphten-liouv` formalism without approximation. Relaxation integration and zero tolerances are both `1e-5`.
- After creating the spin system and basis, the script constructs the relaxation superoperator `R` and reports its matrix elements for the two states.

## Implementation structure

- Specify the magnetic field, isotopes, Zeeman matrices, and coordinates.
- Set relaxation parameters, basis, and accuracy tolerances.
- Create the spin system, construct `R`, and report its action on longitudinal magnetization and the singlet state.
