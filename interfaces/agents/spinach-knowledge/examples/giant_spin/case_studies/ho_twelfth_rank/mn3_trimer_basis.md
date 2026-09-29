# examples/giant_spin/case_studies/ho_twelfth_rank/mn3_trimer_basis.m

- MATLAB implementation: [examples/giant_spin/case_studies/ho_twelfth_rank/mn3_trimer_basis.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/giant_spin/case_studies/ho_twelfth_rank/mn3_trimer_basis.m)

Signature: `[P,msz]=mn3_trimer_basis(spin_system,nstates)`

This helper constructs the reduced spaces used for the (CH6N3)2MnCl4 trimer study. The model has three E6 manganese centres (S=5/2), so the full Zeeman-Hilbert space has 216 states. The returned spaces are built from the isotropic nearest-neighbour exchange Hamiltonian and total S_z, rather than from the full exchange-plus-zero-field-splitting Hamiltonian.

The input spin system must use the zeeman-hilb formalism and contain exactly three E6 spins. Choose `nstates=16` or `nstates=26`, corresponding to Tables 1 and 2 in Section 3.2 of https://arxiv.org/abs/2609.16352. The routine obtains the exchange Hamiltonian by subtracting the unit-field Zeeman operator from the lab-frame Hamiltonian, then adds a 1 mT Zeeman perturbation to resolve the exchange states by total S_z. Energies are referenced to the lowest state and converted to cm^-1 for selection.

For the 16-state space it retains the lowest-energy state in each of the sixteen total-S_z sectors, from -15/2 through +15/2 in unit steps. For the 26-state space it retains all eigenstates with energy strictly below 25 cm^-1 above the ground state. The outputs are `P`, a 216-by-nstates matrix whose columns are the selected states, and `msz`, their total S_z projections. The helper checks the formalism, isotope list, and requested dimension; it does not return the projected Hamiltonian or include zero-field splitting in the basis-selection Hamiltonian.

The resulting spaces are used by the companion level-diagram and pulsed-magnetisation examples. The 16- and 26-state choices and their selection rules follow the cited paper.

Tables 1–2 of the cited paper list multiplet energies above the ground state (cm^-1) at 0 (S=5/2), 12.10 (S=3/2), 16.94 (S=7/2), 24.20, 38.72 (S=9/2), 65.34 (S=11/2), 96.80 (S=13/2), and 133.10 (S=15/2).