# examples/giant_spin/case_studies/mn3_trimer_levels.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/giant_spin/case_studies/mn3_trimer_levels.m`
- Signature: `mn3_trimer_levels()`
- Total lines: 59

## Purpose

Zeeman energy level diagram of the (CH6N3)2MnCl4 molecular crystal, a linear trimer of three S=5/2 manganese ions with isotropic exchange between neighbours and an axial plus rhombic zero-field splitting on every ion, from zero to 10 Tesla in the full 216-state Hilbert space. The lowest levels are the ones that the 16-state and 26-state effective bases of the paper are built to reproduce. Reproduces Fig 5 of https://arxiv.org/abs/2609.16352. Calculation time: seconds.

## Physical / mathematical content

- The same spin system as `mn3_trimer_magn.m`: three `E6` spins with g=2, nearest-neighbour exchange J=-2.42 cm^-1 in the paper's H=-2*J*S1*S2 convention, and D=0.167 cm^-1, E=0.040 cm^-1 on every ion through `zfs2mat`.
- The level crossings and avoided crossings of the lowest thirty levels are where the pulsed-field magnetisation of the companion script changes; the paper's effective bases are chosen to reproduce exactly these levels.

## Numerical / algorithmic content

- `zeeman-hilb` with no basis reduction; the field-free Hamiltonian is the labframe Hamiltonian minus the unit-field Zeeman operator (`sys.magnet=1`), and the Zeeman operator per Tesla is added back scaled by each of 201 fields from 0 to 10 T.
- Eigenvalues of the symmetrised full Hamiltonian at each field are converted from rad/s to cm^-1 with `hz2icm(E/(2*pi))` and plotted relative to the field-free ground state.
