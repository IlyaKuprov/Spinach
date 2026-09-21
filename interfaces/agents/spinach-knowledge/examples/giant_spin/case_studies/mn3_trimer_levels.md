# examples/giant_spin/case_studies/mn3_trimer_levels.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/giant_spin/case_studies/mn3_trimer_levels.m`
- Signature: `mn3_trimer_levels()`
- Total lines: 85

## Purpose

Zeeman energy level diagram of the (CH6N3)2MnCl4 molecular crystal, a linear trimer of three S=5/2 manganese ions with isotropic exchange between neighbours and an axial plus rhombic zero-field splitting on every ion, from zero to 10 Tesla in the full 216-state Hilbert space. The lowest levels are the ones that the 16-state and 26-state effective bases of the paper are built to reproduce. Reproduces Fig 5 of https://arxiv.org/abs/2609.16352 with the axis limits, colours, and the slight vertical offsets that the paper applies to show the overlapping lines; the 16-state and 26-state effective bases of the paper (common eigenstates of the isotropic exchange and the total S_z) come from mn3_trimer_basis.m. Calculation time: seconds.

## Physical / mathematical content

- The same spin system as `mn3_trimer_magn.m`: three `E6` spins with g=2, nearest-neighbour exchange J=-2.42 cm^-1 in the paper's H=-2*J*S1*S2 convention, and D=0.167 cm^-1, E=0.040 cm^-1 on every ion through `zfs2mat`.
- The effective bases reproduce the lowest levels well; where they miss levels of the full space (the paper's regions R2 to R4) the projected diagram departs from the grey full-space lines, which is why the paper's 16-state and 26-state magnetisation curves differ near 10 T.

## Numerical / algorithmic content

- `zeeman-hilb` with no basis reduction; the field-free Hamiltonian is the labframe Hamiltonian minus the unit-field Zeeman operator (`sys.magnet=1`), and the Zeeman operator per Tesla is added back scaled by each of 201 fields from 0 to 10 T; the same Hamiltonian is projected onto the 16-state and 26-state bases before diagonalisation.
- Eigenvalues are converted from rad/s to cm^-1 with `hz2icm(E/(2*pi))`, plotted relative to the field-free ground state of the full space, in grey (all 216 states), black (16 states, raised by 0.15 cm^-1), and red (26 states, raised by 0.3 cm^-1), on the paper's axes of 0 to 10 T and -10 to 10 cm^-1.
