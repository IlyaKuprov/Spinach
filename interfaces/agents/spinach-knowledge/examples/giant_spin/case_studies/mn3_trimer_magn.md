# examples/giant_spin/case_studies/mn3_trimer_magn.m

- Signature: `mn3_trimer_magn()`

## Purpose

Pulsed-field magnetisation of the (CH6N3)2MnCl4 molecular crystal, a linear trimer of three S=5/2 manganese ions with isotropic exchange between neighbours and an axial plus rhombic zero-field splitting on every ion, at 0.6 K under a 50 T/ms sweep to 10 T with spin-phonon relaxation in the generalised Lindblad form of Saito and Miyashita. As in the paper, the dynamics runs in the 16-state and 26-state effective bases of mn3_trimer_basis.m, with the Hamiltonian, the Zeeman operator, and the observable projected onto each basis; the thermal equilibrium magnetisation of the full 216-state space is plotted for comparison. Reproduces Fig 6 of https://arxiv.org/abs/2609.16352 with the colours and axis limits of that paper. Calculation time: minutes.

## Physical / mathematical content

- Three `E6` spins with g=2, nearest-neighbour isotropic exchange J=-2.42 cm^-1 in the paper's H=-2*J*S1*S2 convention (`-2*icm2hz(-2.42)*eye(3)` in Spinach's S1*A*S2 form), and `zfs2mat(icm2hz(0.167),icm2hz(0.040),0,0,0)` on every ion; the sweep runs in the frame of the zero-field splitting tensors.
- Both effective bases give a plateau near 3.35 mu_B below the equilibrium value near 5 mu_B, a sudden rise near 4.9 T (the paper's sizable sudden increase at about 5 T) to a second plateau near 3.9 mu_B, and a rise beyond 9 T; the two bases part near 10 T, where the 26-state basis lacks the |S_z|>7/2 states and the 16-state basis lacks the thermally relevant low states, as the paper discusses.
- The spin-phonon coupling operator has unit elements between basis states whose total S_z differs by one, as in qdmag; super-Ohmic bath (`phonon_alpha=2`) with lambda=10 cm^-1 and I0=1e-14 ps/rad; the observable is the total moment -2*S_z in Bohr magnetons.

## Numerical / algorithmic content

- `zeeman-hilb`; the full Hamiltonian at 1 T (`sys.magnet=1`) and the Zeeman operator per Tesla are projected with `P'*H*P` onto each basis from `mn3_trimer_basis`, symmetrised, and handed to `pulsed_field` directly with `parameters.hzeeman` set, bypassing the `crystal` context; the spin system carries the `labframe` assumption through `assume` because `pulsed_field` checks it. 20 ns stairs, 10^4 stairs, one record every 100 stairs.
- The equilibrium curve at every recorded field is `equilibrium(spin_system,H)` on the full-space Hilbert space Hamiltonian H0+B*Z, read out with `real(hdot(coil,rho))`; plotted in blue with the 16-state sweep in red and the 26-state sweep in green on the paper's axes of 0 to 10 T and 0 to 6 mu_B.
