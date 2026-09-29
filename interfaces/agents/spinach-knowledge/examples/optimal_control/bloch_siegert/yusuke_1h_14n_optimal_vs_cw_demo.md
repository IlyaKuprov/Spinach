# examples/optimal_control/bloch_siegert/yusuke_1h_14n_optimal_vs_cw_demo.m

> **Historical example — no longer shipped in current Spinach.** This page describes a deleted example only; do not try to run it or assume its source path exists in the current package.

- Historical function signature: yusuke_1h_14n_optimal_vs_cw_demo()

## Historical model and setup

The script set up a reduced two-spin model: one observed 1H and one controlled 14N. Its header describes the field as corresponding to 800 MHz 1H; the code sets sys.magnet=18.8 and a scalar-coupling entry of 1500. It treats the coupling as an effective reduced-model surrogate, not a literal solid-state Hamiltonian. Controls act on 14N, while the target is preservation of the 1H transverse operators.

The script seeded its random-number generator with 1. The phase-only control used a nominal 14N RF value of 20 kHz and 120 pulse elements of 10 μs each (1.2 ms total from those settings). Optimisation used a square-fidelity objective and L-BFGS, with a maximum of 40 iterations and Bloch–Siegert correction enabled on the driven 14N channel. The training ensemble used seven offsets from −12 to +12 kHz and B1 scales 0.95, 1.00, and 1.05. The script initialised the phase sequence with a repeating four-phase pattern before calling its GRAPE phase optimiser.

It compared the optimised phase waveform with a zero-phase, constant-amplitude CW-like baseline of the same duration and nominal RF value. Evaluation sampled 61 14N offsets from −20 to +20 kHz and nine B1 scales from 0.90 to 1.10, then formed offset/B1 profiles and training-ensemble scores. The script is written to plot these evaluations and print mean and minimum training fidelities for each waveform; no numerical scores or successful-run claim are asserted here.

The source describes the example as inspired by offset-tolerant 14N decoupling work by Nehra, Agarwal, and Nishiyama, but gives no complete bibliographic reference or DOI; none is inferred here.
