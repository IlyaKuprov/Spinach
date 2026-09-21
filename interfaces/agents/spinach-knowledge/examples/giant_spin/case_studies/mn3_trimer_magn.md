# examples/giant_spin/case_studies/mn3_trimer_magn.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/giant_spin/case_studies/mn3_trimer_magn.m`
- Signature: `mn3_trimer_magn()`
- Total lines: 86

## Purpose

Pulsed-field magnetisation of the (CH6N3)2MnCl4 molecular crystal, a linear trimer of three S=5/2 manganese ions with isotropic exchange between neighbours and an axial plus rhombic zero-field splitting on every ion, at 0.6 K under a 50 T/ms sweep to 10 T with spin-phonon relaxation in the generalised Lindblad form of Saito and Miyashita. The full 216-state Hilbert space is used; the paper solves the same problem in 16-state and 26-state effective bases. The thermal equilibrium magnetisation is plotted for comparison. Reproduces Fig 6 of https://arxiv.org/abs/2609.16352. Calculation time: minutes.

## Physical / mathematical content

- Three `E6` spins with g=2, nearest-neighbour isotropic exchange J=-2.42 cm^-1 in the paper's H=-2*J*S1*S2 convention (`-2*icm2hz(-2.42)*eye(3)` in Spinach's S1*A*S2 form), and `zfs2mat(icm2hz(0.167),icm2hz(0.040),0,0,0)` on every ion; the sweep runs in the frame of the zero-field splitting tensors.
- Below about 4.8 T the magnetisation sits on a plateau near 3.35 mu_B, well under the equilibrium value near 4.9 mu_B; near 4.9 T it rises suddenly, the feature the paper describes as a sizable sudden increase at about 5 T, and it reaches about 5.6 mu_B at 10 T against 5.0 at equilibrium, because the phonon bath at 0.6 K cannot keep up with a 50 T/ms sweep.
- The spin-phonon coupling operator has unit elements between product states whose total S_z differs by one; super-Ohmic bath (`phonon_alpha=2`) with lambda=10 cm^-1 and I0=1e-14 ps/rad; the observable is the total moment -2*S_z in Bohr magnetons.

## Numerical / algorithmic content

- `zeeman-hilb` with no basis reduction (216 states), `crystal` context with the `labframe` assumption and `parameters.needs={'zeeman_op'}`, `sys.magnet=1`; 20 ns stairs, 10^4 stairs, one record every 100 stairs. Total S_z comes from `operator(spin_system,'Lz','E6')`.
- The equilibrium curve at every recorded field is `equilibrium(spin_system,H)` on the Hilbert space Hamiltonian H0+B*Z, where H0 is the labframe Hamiltonian minus the unit-field Zeeman operator and Z is the Zeeman operator per Tesla; the expectation value is `real(hdot(coil,rho))`.
