# Cavity-assisted nuclear pumping: native Liouville-space example

Run `result=cavity_dnp_minimal(3)` and `validation=validate_cavity_dnp()` after adding Spinach and this directory to the MATLAB path. Neither function writes files or opens figures.

## Physical scope

This is a **four-level electron–spin-1/2 nucleus analogue** of the shelving–nuclear-swap–radiative-reset process in Fig. 1 of Lorriaux et al., arXiv:2610.19180v1. It is not a simulation of their low-field Bi:Si donor, whose strongly mixed electron–nuclear eigenstates and 20-level structure are not retained. It does not predict that sample's population, nuclear polarization, timescale, or collective enhancement.

The secular hyperfine coupling is 10 kHz, the cavity detuning is 5 kHz, the exchange coupling is 100 Hz, and the photon linewidth is 2 kHz. These illustrative parameters put the nuclear-alpha electron transition on resonance and give well-resolved cavity-mediated reset. They are not the paper's fitted parameters. The common rotating frame and rotating-wave approximation are requested with `assume(spin_system,'cavity')`; the electron–nuclear coupling strength is then set to `zz` to retain only SzIz. The cavity bath is at zero temperature; intrinsic electron/nuclear relaxation, ensemble disorder, and collective effects are absent.

Initially the electron is in its lower Zeeman level, the nucleus is unpolarized, and the cavity is in vacuum. An ideal nuclear-alpha-selective electron pi pulse shelves half the population. An electron-ground-state-selective nuclear pi pulse moves the other half into the target nuclear state. Native Jaynes–Cummings exchange and photon loss then reset the shelved electron into that same target state. The target population approaches one from one half, and nuclear Iz approaches +1/2 from zero. Ideal pulses neglect drift and relaxation during the pulse and do not model the experimental WURST waveforms.

## Native Spinach implementation

There is exactly one formalism: `sphten-liouv`, with the complete basis. Native `paranoia` precision settings prevent propagator chopping from contaminating Fock-cutoff comparisons. The particles are `E`, `1H`, and `Cn`. The entire Hamiltonian comes from `inter.coupling.scalar`, `inter.modes.frqs`, and `inter.modes.exchange`; the dissipator comes from `inter.modes.linewidths`, `inter.temperature`, and `relaxation`.

Initial states and observables use native `state`/`operator` syntax with `ZL1`, `ZL2`, `BL1`, `BLn`, `Lz`, and `N`. Observable dual vectors are `operator(...,'left')'*unit`, giving physical Tr(O rho) without assuming a Euclidean coordinate metric. Spinach numbers spin levels from the bottom: `ZL1` is the spin-1/2 lower level, while `BL1` is bosonic vacuum. Initial trace normalization is obtained with the native identity state. Conditional pulse generators use native `operator` product syntax; `step` applies each ideal pulse and `evolution` computes the reset trajectory. There are no handwritten matrices, dissipators, propagators, tensor products, spectral projectors, or formalism conversions.

## Validation

The validator compares cavity dimensions 3, 4, and 5 and halves the time step. It checks trace, pulse population swaps, nonnegative reported populations, negligible top-level occupation, a zero-exchange control, lossless excitation conservation and resonant Rabi oscillations, and the exact resolved damped Jaynes–Cummings photon transient. Analytical curves are validation references only, never substituted into the simulated dynamics.
