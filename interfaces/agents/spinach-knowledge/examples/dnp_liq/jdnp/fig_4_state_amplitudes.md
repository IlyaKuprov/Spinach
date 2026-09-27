# examples/dnp_liq/jdnp/fig_4_state_amplitudes.m

- Signature: `fig_4_state_amplitudes()`
- Reference: [Concilio et al., *Phys. Chem. Chem. Phys.* (2022)](https://doi.org/10.1039/d1cp04186j)
- Calculation time: seconds

## Purpose

Propagates a liquid-state radical-pair model under microwave irradiation and plots the populations of the electron-pair singlet and triplet states resolved by nuclear-spin projection. The example illustrates how the singlet-alpha and singlet-beta populations evolve differently, alongside the triplet populations and nuclear magnetisation.

## Model and setup

The script loads the system from `system_specification()`, sets a 14.08 T field, and drives the electron spins with a 250 kHz microwave field. The microwave offset is set from the trityl and free-electron resonance frequencies. It changes the electron-proton scalar coupling using the electron and proton Zeeman frequencies, then constructs the spin system in the supplied basis.

## Propagation and output

The Hamiltonian is built for ESR conditions, combined with the relaxation superoperator and microwave terms, and propagated from the isotropic thermal-equilibrium state. The `evolution` call requests 700 steps at 1 ms spacing in multichannel mode for the state operators assembled in the script. The resulting figure has three panels: triplet populations for alpha and beta nuclear projections, the two singlet populations, and the nuclear Lz signal. Curves are plotted as real parts against time.
