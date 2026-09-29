# examples/singlet_states/carbon_singlet.m

- Signature: `carbon_singlet()`

## Purpose

Construct a Redfield relaxation superoperator for the two triple-bond carbons in cis-dimethylbut-2-ynedioate, then report its matrix elements for normalised longitudinal magnetisation and the two-spin singlet state. The source says the magnetic parameters were computed with DFT.

## Spin pair and relaxation model

The system is a pair of 13C spins at 14.1 T. The source supplies Zeeman matrices `[[29.13, 0, 0], [0, 253.00, -5.64], [0, -35.61, 38.78]]` and `[[29.13, 0, 0], [0, 253.00, 5.64], [0, 35.61, 38.78]]`, with coordinates (0, 0.609, 0.298) and (0, -0.609, 0.298). It selects Redfield relaxation, zero equilibrium, lab-frame retention, and a 100 ps correlation time. The basis is the full `sphten-liouv` basis (`approximation='none'`); the relaxation integration and zero tolerances are both `1e-5`.

## Observable and limits

After constructing `R`, the example normalises `Sz` and the singlet operator and reports `Sz'*R*Sz` and `S'*R*S`. These are the reported relaxation-superoperator matrix elements; the script does not evolve a prepared state or report numerical lifetimes. It contains no storage interval, imaging or gradient model.

The source comment gives a calculation time of seconds.

Source: [examples/singlet_states/carbon_singlet.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/singlet_states/carbon_singlet.m)
