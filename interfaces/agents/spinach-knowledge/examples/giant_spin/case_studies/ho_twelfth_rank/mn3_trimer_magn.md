# examples/giant_spin/case_studies/ho_twelfth_rank/mn3_trimer_magn.m

- MATLAB implementation: [examples/giant_spin/case_studies/ho_twelfth_rank/mn3_trimer_magn.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/giant_spin/case_studies/ho_twelfth_rank/mn3_trimer_magn.m)

Signature: `mn3_trimer_magn()`

This example calculates pulsed-field magnetisation for the (CH6N3)2MnCl4 linear trimer: three S=5/2 E6 spins, g=2, nearest-neighbour isotropic exchange J=-2.42 cm^-1 in the paper's H=-2*J*S1.S2 convention, and axial/rhombic zero-field splitting D=0.167 cm^-1 and E=0.040 cm^-1 on every ion. The calculation uses the same tensor orientation and effective bases as the companion level diagram and targets the Fig. 6 comparison in https://arxiv.org/abs/2609.16352.

The spin system is initialised with a 1 T reference magnet. At 0.6 K, a super-Ohmic spin-phonon bath is used in the generalised Lindblad treatment of Saito and Miyashita (alpha=2; lambda=10 cm^-1 and I0=1e-14 ps/rad as specified in the paper). The sweep is 50 T/ms, from 0 to 10 T: a 20 ns step for 10,000 steps, with output every 100 steps. For each of the 16- and 26-state bases, the Hamiltonian, Zeeman operator per Tesla, and moment operator are projected into the basis before the pulsed-field simulation. The phonon coupling connects basis states whose total S_z labels differ by one. The measured quantity is the total moment along Z, -2*S_z, reported in Bohr magnetons (mu_B).

The figure compares the two reduced-space dynamical curves with the full 216-state thermal-equilibrium magnetisation. It uses a 0-10 T field range and a 0-6 mu_B vertical range; equilibrium, 16-state, and 26-state traces are blue, red, and green. The 16-state and 26-state calculations use different reduced bases from `mn3_trimer_basis`; they should not be confused with propagation in the full 216-state space.
