# examples/giant_spin/case_studies/dimer_exchange_types.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/giant_spin/case_studies/dimer_exchange_types.m`
- Signature: `dimer_exchange_types()`
- Total lines: 92

## Purpose

Pulsed-field magnetisation of a dimer of two S=1/2 spins with four types of exchange coupling tensor: isotropic, two anisotropic, and antisymmetric, at 0.2 K under a 10 T/ms sweep to 1 T, with spin-phonon relaxation in the generalised Lindblad form of Saito and Miyashita. The out-of-equilibrium curves are compared with the thermal equilibrium magnetisation. Reproduces Fig 4 of https://arxiv.org/abs/2609.16352. Calculation time: minutes.

## Physical / mathematical content

- Two electrons with g=2 and one exchange tensor at a time: isotropic 0.2 cm^-1, J_xx only, J_zz only, and a purely antisymmetric tensor. The paper writes H=-2*S1*J*S2, Spinach writes S1*A*S2 with A in hertz, so `inter.coupling.matrix{1,2}=-2*icm2hz(J)`.
- The isotropic and J_zz tensors commute with the Zeeman term, so the sweep cannot move population between Zeeman levels and the magnetisation stays at zero while the equilibrium curve grows; J_xx and the antisymmetric tensor mix the levels and the spin-phonon dissipator drives the magnetisation towards, but behind, equilibrium (J_xx plateaus near 0.61 mu_B against 1.98 at equilibrium, antisymmetric reaches 1.61 against 1.90 at 1 T).
- The spin-phonon coupling operator has unit elements between product states whose total S_z differs by one; the bath is super-Ohmic (`phonon_alpha=2`) with the paper's lambda^2*I0 (lambda=10 cm^-1, I0=1e-10 ps/rad) converted to rad/s units; the observable is the total moment -2*S_z in Bohr magnetons.

## Numerical / algorithmic content

- `zeeman-hilb` formalism, `crystal` context with the `labframe` assumption and `parameters.needs={'zeeman_op'}`, `sys.magnet=1` so that `pulsed_field` receives the Zeeman operator per Tesla; 10 ns stairs, 10^4 stairs, one record every 10 stairs. The spin system is re-created for every tensor, and the total S_z (`operator(spin_system,'Lz','E')`), the coupling mask, and the coil are built after `basis` inside the loop.
- The equilibrium curve at every recorded field is `equilibrium(spin_system,H)` on the Hilbert space Hamiltonian H0+B*Z, where H0 is the labframe Hamiltonian minus the unit-field Zeeman operator and Z is the Zeeman operator per Tesla from `hamiltonian(assume(spin_system,'labframe','zeeman'))`; the expectation value is `real(hdot(coil,rho))`. The equilibrium values are stored in `answers{n}.obs_eq`.
