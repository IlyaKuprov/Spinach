# examples/giant_spin/case_studies/mn3_trimer_basis.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/giant_spin/case_studies/mn3_trimer_basis.m`
- Signature: `[P,msz]=mn3_trimer_basis(spin_system,nstates)`
- Total lines: 72

## Purpose

Effective basis of the (CH6N3)2MnCl4 trimer case studies, built as in Section 3.2 of https://arxiv.org/abs/2609.16352: the common eigenstates of the isotropic exchange Hamiltonian and the total S_z, obtained by diagonalising the exchange Hamiltonian with a perturbative 1 mT Zeeman term, are selected either as the lowest state of every total S_z projection (16 states, Table 1 of the paper) or as every state within 25 cm^-1 of the ground state (26 states, Table 2 of the paper).

```matlab
[P,msz]=mn3_trimer_basis(spin_system,nstates)
```

## Parameters / inputs

- spin_system -Spinach spin system of the trimer, three E6 spins in zeeman-hilb formalism
- nstates -16 or 26, the basis set of the paper to return

## Outputs

- P -matrix with the basis states in columns, 216 x nstates
- msz -column of total S_z projections of the basis states

## Physical / mathematical content

- The isotropic exchange Hamiltonian is the isotropic part returned by `hamiltonian` minus the unit-field Zeeman operator; adding a 1 mT Zeeman term lifts the S_z degeneracy of each exchange multiplet without mixing multiplets, so `eig` returns common eigenvectors of exchange and total S_z, whose projections are read off as rounded diagonal elements of the transformed S_z.
- With the paper's parameters the multiplets sit at 0 (S=5/2), 12.10 (S=3/2), 16.94 (S=7/2), 24.20, 38.72 (S=9/2), 65.34 (S=11/2), 96.80 (S=13/2), and 133.10 cm^-1 (S=15/2) above the ground state, which reproduces Tables 1 and 2 of the paper: the 16-state basis takes the lowest state of every S_z from -15/2 to 15/2, the 26-state basis every state below 25 cm^-1.

## Numerical / algorithmic content

- `mn3_trimer_levels.m` diagonalises the full Hamiltonian projected onto each basis, and `mn3_trimer_magn.m` runs `pulsed_field` directly with the projected Hamiltonian, Zeeman operator, and observable, and with the coupling operator built from the S_z labels of the basis states, as qdmag does.
