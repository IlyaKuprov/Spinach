# examples/esr_sol_pulsed/holeburn_nitroxide_powder.m

- MATLAB implementation: [examples/esr_sol_pulsed/holeburn_nitroxide_powder.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/esr_sol_pulsed/holeburn_nitroxide_powder.m)

[MATLAB example](../../../../../examples/esr_sol_pulsed/holeburn_nitroxide_powder.m) · [holeburn sequence helper](../../../../../experiments/holeburn.m)

- Signature: `holeburn_nitroxide_powder()`

## Aim and spin system

This example compares three powder-averaged nitroxide hole-burning conditions: a frequency-swept soft pulse, a fixed-frequency soft pulse, and a zero-power reference. The two spins are an electron (`E`) and `14N` at 3.5 T. The electron g tensor is diagonal, `diag(2.01045, 2.00641, 2.00211)`. The electron–nitrogen coupling tensor is supplied as `1e7 * [1.2356 0 0.6322; 0 1.1266 0; 0.6322 0 8.2230]`; the source does not annotate the matrix unit. The exact spherical-tensor Liouville basis is used (`approximation='none'`); trajectory-level SSR algorithms are disabled.

## Pulse and acquisition protocol

The `holeburn` helper models the soft pulse by Fokker–Planck propagation and follows it with an ideal hard `pi/2` observation pulse before acquisition. The shaped pulse uses rank 2, phase `pi/2` rad, method `expv`, 100 pulse steps of `1e-9` s each, and power `2*pi*10e6` rad/s at each step. In the chirp condition (A), the 100 carrier-frequency entries run linearly from `-350e6` to `-250e6` Hz. In the fixed-frequency condition (B), all 100 entries are `-300e6` Hz. The reference (C) has zero power and zero carrier-frequency entries. Thus the pulse power and step durations are shared by A and B; the swept versus fixed carrier is the comparison.

The initial state is `Lz` and the coil operator is `L+` for the electron; no spins are decoupled. Acquisition uses offset `-2e8` (unit not annotated), sweep width `8e8` Hz, 64 points, zero-fill to 512, the `rep_2ang_3200pts_sph` powder grid, no derivative, and no axis inversion; the display unit is MHz. The helper documentation specifies Hz for pulse frequency, rad/s for pulse power, and seconds for pulse duration; `acquire` specifies sweep width in Hz.

## Observable and output

Each FID is exponentially apodised with parameter 6 and Fourier transformed. The real spectra are overlaid in a figure: red for the chirp, blue for the fixed-frequency soft pulse, and black for the reference. The source does not save a data or figure file and estimates a calculation time of seconds. This is a simulated powder-averaged signal with an ideal observation pulse; it is not a measured spectrum.
