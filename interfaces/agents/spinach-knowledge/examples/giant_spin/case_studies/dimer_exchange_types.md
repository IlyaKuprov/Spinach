# examples/giant_spin/case_studies/dimer_exchange_types.m

- Signature: `dimer_exchange_types()`
- Current source: `examples/giant_spin/case_studies/ho_twelfth_rank/dimer_exchange_types.m`. The H1 path above is the knowledge-base entry key, not the current source location.

## Purpose

This example simulates pulsed-field magnetisation for a dimer of two S=1/2 spins with four exchange tensors: isotropic, two anisotropic, and antisymmetric. At 0.2 K, the field sweeps to 1 T at 10 T/ms. Spin-phonon relaxation is treated with the generalised Lindblad form of Saito and Miyashita, and the resulting non-equilibrium curves are compared with thermal-equilibrium magnetisation. The calculation reproduces Fig 4 of https://arxiv.org/abs/2609.16352, including its colours, line styles, and axis limits. Calculation time: minutes.

## Physical / mathematical content

- The two electron spins have g=2. One exchange tensor is used at a time: isotropic 0.2 cm^-1, J_xx only, J_zz only, or purely antisymmetric. The paper uses H=-2*S1*J*S2, whereas Spinach uses S1*A*S2 with A in hertz; accordingly, `inter.coupling.matrix{1,2}=-2*icm2hz(J)`.
- The isotropic and J_zz tensors commute with the Zeeman term, so they do not coherently mix Zeeman levels. Population transfer then comes only from the spin-phonon dissipator, which couples product states whose total S_z differs by one. At the paper's spectral density, magnetisation remains near 0.006 and 0.003 mu_B, respectively, at 1 T, compared with 2.0 at equilibrium. J_xx and the antisymmetric tensor mix the levels; their magnetisation follows equilibrium with a lag. J_xx plateaus near 0.61 mu_B versus 1.98 at equilibrium, while the antisymmetric case reaches 1.61 versus 1.90 at 1 T.
- The spin-phonon coupling operator has unit elements between product states whose total S_z differs by one. The bath is super-Ohmic (`phonon_alpha=2`) and uses the paper's lambda^2*I0, with lambda=10 cm^-1 and I0=1e-10 ps/rad, converted to rad/s units. The observable is the total moment -2*S_z in Bohr magnetons.

## Numerical / algorithmic content

- The calculation uses the `zeeman-hilb` formalism and a `crystal` context with the `labframe` assumption and `parameters.needs={'zeeman_op'}`. Setting `sys.magnet=1` supplies `pulsed_field` with the Zeeman operator per Tesla. The sweep uses 10 ns steps for 10^4 steps, with output every 10 steps. The spin system is recreated for each tensor; after `basis` is set, the loop constructs total S_z (`operator(spin_system,'Lz','E')`), the coupling mask, and the coil.
- At each recorded field, the equilibrium curve is computed from `equilibrium(spin_system,H)` on the Hilbert-space Hamiltonian H0+B*Z. Here H0 is the lab-frame Hamiltonian minus the unit-field Zeeman operator, and Z is the Zeeman operator per Tesla from `hamiltonian(assume(spin_system,'labframe','zeeman'))`. The expectation value is `real(hdot(coil,rho))`; equilibrium values are stored in `answers{n}.obs_eq`.
- A single panel follows the paper's styling: each tensor's equilibrium curve is solid and its sweep curve dashed, in red, blue, green, and orange for tensors 1 to 4. The axes span 0 to 1 T and 0 to 2 mu_B. A transparent legend on the right lists "J_n^dimer Equilibrium" and "J_n^dimer QME".
