# examples/quantum_tech/spin_cavity_purcell_effect.m

- Signature: `spin_cavity_purcell_effect()`
- Source: [`examples/quantum_tech/spin_cavity_purcell_effect.m`](https://github.com/IlyaKuprov/Spinach/blob/main/examples/quantum_tech/spin_cavity_purcell_effect.m)

## Purpose and physical model

A Liouville-space model of cavity-induced electron-spin relaxation in the EPR Purcell regime. It combines coherent Jaynes–Cummings exchange with cavity damping and compares relaxation rates across spin–cavity detuning and cavity loss. The declared isotopes are `{'E','C3'}`: one generic Spinach electron spin and a three-level cavity mode; the source sets the magnetic field and temperature to zero. No particular chemical defect or defect-specific nuclear spin is identified.

The spin–cavity exchange parameter is `0.35e6` Hz. The code samples angular detuning from `-2*pi*8e6` to `+2*pi*8e6` rad/s at `301` points. Its four cavity loss settings are `kappa=2*pi*[2e6 4e6 8e6 16e6]` rad/s, labelled in the plots as `kappa/(2*pi)=2, 4, 8, 16` MHz. The helper builds the resonant cavity at zero frequency, converts the linewidth to the Spinach Hz convention, sets the exchange coupling, and constructs the damping superoperator with `relaxation`. The basis is `zeeman-liouv` with no approximation.

There is no pulse sequence, independent coherent-drive term, or simulated defect spectrum: the model isolates an electron spin transition using its `Lz` superoperator and sweeps its detuning from the damped cavity. Cavity loss enters through the relaxation superoperator. This is an EPR-motivated generic spin–cavity model, not a calculation of a measured defect line.

## Rate extraction and plotted observables

For each detuning and loss setting, the source forms the dissipative Liouvillian `L=-1i*(Hjc+detuning*spin_ham)+R`. It takes the slow retained eigenmode's decay rate and multiplies its amplitude-decay rate by two to report a population-relaxation rate. At resonance and the reference loss setting it also checks against the exact two-level expression `kappa/2-sqrt(kappa^2/4-4*g_rad^2)`, where `g_rad=2*pi*coupling`; these are source assertions, not independently rerun results here.

The left panel plots that calculated rate versus detuning, displayed in MHz and kHz, for the four loss settings. The right panel propagates an initially excited electron spin with the cavity in `BL1` and plots the spin-excitation survival for detunings `2*pi*[0 2e6 6e6]` rad/s at the reference `kappa/(2*pi)=4` MHz, over 0–40 μs with `250` time samples. These are model trajectories and Liouvillian eigenvalue calculations, not measurements of a physical cavity.
