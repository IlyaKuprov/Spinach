# examples/quantum_tech/circuit_qed/cavity_dual_rail_safe.m

- Signature: `cavity_dual_rail_safe()`

## Purpose

Calculates flux-noise dephasing rates for two dual-rail qubits whose transmon-coupled rails share one flux-tunable transmon, and tests their suppression by a SAFE drive. The example follows Sec. 4.4.2 and Fig. 4.4(c) of Yunwei Lu's 2026 Northwestern University PhD thesis. The driven rates are computed from eigenvalues of the rotating-frame Hamiltonian in Eq. (4.16) at fixed drive amplitude as functions of transmon-drive detuning; the SAFE relations are Eqs. (4.93) and (4.95).

## Physical / mathematical content

The Hilbert-space model contains a three-level transmon and two three-level cavity rails. Their exchange couplings are 50 and 60 MHz, detunings are 1.414 and 2.0 GHz, and the transmon anharmonicity is -200 MHz. The transmon drive amplitude is 10 MHz; the code scans drive detunings from -20 to -50 MHz in 0.5 MHz steps. Driven rail dephasing rates are calculated from flux sensitivities of the dressed single-photon transition frequencies; undriven rates use the bare dispersive sensitivities.

## Numerical / algorithmic content

The flux sensitivity is estimated by a centred 1 kHz frequency difference. For each rail, the script locates the minimum rate and chooses a common operating point from the mean of the two minimum indices. It checks that both minima lie inside the scan window, are within 5 MHz of each other, and that the common point suppresses each undriven rate by at least a factor of five. The noise amplitude is `1e-5` flux quanta and the transmon frequency sensitivity is `2*pi*6e9` rad/s per flux quantum.

## Implementation structure

The script builds the transmon-plus-two-rail Hamiltonian and flux-noise operator, evaluates dressed eigenvalues over the detuning scan, computes driven and bare rates, validates the common sweet-spot window and suppression factors, and plots both rails' rates together with the undriven reference lines.
