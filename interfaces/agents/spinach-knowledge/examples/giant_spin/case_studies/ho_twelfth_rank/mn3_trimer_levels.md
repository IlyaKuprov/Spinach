# examples/giant_spin/case_studies/ho_twelfth_rank/mn3_trimer_levels.m

- MATLAB implementation: [examples/giant_spin/case_studies/ho_twelfth_rank/mn3_trimer_levels.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/giant_spin/case_studies/ho_twelfth_rank/mn3_trimer_levels.m)

Signature: `mn3_trimer_levels()`

This example plots the Zeeman levels of the (CH6N3)2MnCl4 molecular crystal, modelled as three E6 manganese ions with S=5/2. The 216-state model has isotropic exchange between neighbours, axial and rhombic zero-field splitting on each ion, and g=2 on all three spins. It is the calculation behind Fig. 5 of https://arxiv.org/abs/2609.16352.

The exchange parameter is J=-2.42 cm^-1 in the paper's H=-2*J*S1.S2 convention. Each ion has D=0.167 cm^-1 and E=0.040 cm^-1, with the zero-field-splitting tensor orientation set to zero. The spin system is initialised with a 1 T reference magnet; the field-free Hamiltonian is separated from the unit-field Zeeman operator and evaluated from 0 to 10 T at 201 equally spaced fields. The same full Hamiltonian is diagonalised in the full 216-state space and after projection into the 16- and 26-state bases from `mn3_trimer_basis.m`.

The plotted energies are converted to cm^-1 and referenced to the full-space ground-state energy at zero field. The figure shows the full spectrum in grey, the 16-state levels in black, and the 26-state levels in red; the reduced-basis traces are shifted upward by 0.15 and 0.30 cm^-1, respectively, solely to make coincident lines visible. The plot limits are 0-10 T and -10 to 10 cm^-1. These offsets are graphical, not changes to the Hamiltonian or calculated energies. The reduced spaces are intended to reproduce the low levels; the paper discusses departures in regions R2-R4 where levels are absent from the truncated descriptions.