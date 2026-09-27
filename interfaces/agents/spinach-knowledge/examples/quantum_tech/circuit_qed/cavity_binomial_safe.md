# examples/quantum_tech/circuit_qed/cavity_binomial_safe.m

- Signature: `cavity_binomial_safe()`

## Purpose

Evaluates SAFE-drive suppression of flux-noise dephasing and decoherence for a binomial cavity code. The code states are `|0L>=(|0>+|4>)/sqrt(2)` and `|1L>=|2>`. This example follows Sec. 4.4.1 and Fig. 4.4(a,b) of Yunwei Lu's 2026 Northwestern University PhD thesis; it refers to Eqs. (4.93) and (4.94).

## Physical / mathematical content

A three-level flux-tunable transmon is dispersively coupled to a five-level cavity. The model uses transmon anharmonicity `-67 MHz`, coupling `86 MHz`, and detuning `1.414 GHz`, giving the stated 0.5 MHz dispersive shift. The SAFE drive amplitude is 10 MHz and its detuning is scanned from -20 to -80 MHz. Flux-noise dephasing rates for the listed code- and error-space coherences are obtained from flux sensitivities of the dressed cavity transition frequencies; the common operating detuning is selected from the median of the individual rate minima.

## Numerical / algorithmic content

The script uses flux-noise amplitude `1e-5`, a 50 MHz ultraviolet cutoff, transmon and cavity relaxation times of 50 microseconds and 20 milliseconds, and 100 trajectories. It propagates a Lindblad master equation for 30,000 steps of 10 ns (300 microseconds), using propagators tabulated on 201 noise-frequency values and recording infidelity every 300 steps. It compares driven and undriven dephasing rates and decoherence-only infidelity, and computes the final cavity Wigner functions on a 71-by-71 grid. Runtime checks require at least fivefold rate suppression and threefold infidelity reduction.

## Implementation structure

The script constructs dressed-state Hamiltonians for the detuning scan, builds the Liouville-space drift, drive, and noise generators, synthesises the noise trajectories, and averages propagated states. It then traces out the transmon, checks the infidelity criterion, computes and normalises the Wigner functions, and plots the rates, infidelity, and cavity phase-space distributions.
