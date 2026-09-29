# examples/quantum_tech/circuit_qed/cavity_binomial_safe.m

Source: https://github.com/IlyaKuprov/Spinach/blob/main/examples/quantum_tech/circuit_qed/cavity_binomial_safe.m

## Purpose and model

This example studies Stark-assisted flux-noise evasion (SAFE) for the binomial cavity code `|0L> = (|0> + |4>)/sqrt(2)` and `|1L> = |2>`, as discussed in Section 4.4.1 and Figure 4.4(a,b) of Yunwei Lu's 2026 Northwestern University PhD thesis. It computes noise-induced dephasing estimates and compares driven and undriven decoherence in a finite-dimensional transmon–cavity model; it does not report a device measurement.

The transmon is truncated to three levels (`T3`) and the cavity to five (`C5`), with zero-temperature, zero-frame-frequency modes and no basis approximation. The drift contains the transmon self-anharmonicity and a dispersive cross-Kerr number-number shift between transmon and cavity; the off-resonant SAFE drive Stark-dresses the transmon and changes the dressed cavity-transition flux slopes. The flux-noise operator combines transmon number, cavity number, and their product, weighted by the corresponding frequency sensitivities. The transmon anharmonicity is −67 MHz. Its exchange coupling to the cavity is 86 MHz and the detuning is 1.414 GHz; the source states that these give a 0.5 MHz dispersive-shift magnitude. The transmon and cavity relaxation times are 50 μs and 20 ms. The SAFE drive is 10 MHz in amplitude, with transmon-drive detunings scanned from −20 to −80 MHz in 1 MHz steps.

## Flux-noise calculation and dynamics

For the code and error-space coherence pairs (2,0), (4,2), (4,0), (3,0), (4,3), (2,1), and (3,1), the script estimates dephasing from dressed cavity-transition flux sensitivities, using symmetric frequency perturbations of `±2*pi*1e3` rad/s and following thesis Eq. (4.93). It compares the rates with the drive on and off and selects a common operating detuning from the rate minima. The model uses flux amplitude `1e-5` flux quanta, transmon frequency sensitivity `2*pi*6e9` rad/s per flux quantum, a 50 MHz upper noise cutoff, and a 10 ns step. The infrared cutoff follows from a `2^19`-point noise grid and is about 190.7 Hz. The initial logical `|+L>` state has cavity amplitudes `(1, 0, sqrt(2), 0, 1)/2` in the five-level Fock basis. These are simulated noise and model parameters, not measured device values.

The logical `|+L>` state is propagated under a Lindblad generator with transmon/cavity relaxation and 1/f flux-noise trajectories. The comparison uses 100 trajectories, 30,000 steps (300 μs total), and records every 300 steps. A table of propagators is evaluated on 201 flux-frequency offsets. The code compares driven and undriven decoherence-only infidelity (thesis Eq. 4.94) against ideal closed-system evolution, traces out the transmon to obtain the cavity state, and plots cavity Wigner functions on a 71 × 71 phase-space grid. Its runtime assertions require fivefold suppression of every listed rate and a threefold reduction in final infidelity; these are pass/fail thresholds, not numerical results quoted here.

The script implements the thesis model and does not establish experimental protection or convergence beyond its stated truncation, noise bandwidth, time step, and trajectory count.
