# examples/optimal_control/bloch_siegert/yusuke_1h_14n_optimal_vs_cw_demo.m

- Signature: `yusuke_1h_14n_optimal_vs_cw_demo()`

## Status

Historical example: the MATLAB source was removed in commit `b4f03f29`. This note describes the source at its parent revision; the example is not available at the current source path and should not be treated as a runnable current demo.

## Purpose and model

The deleted script is a reduced heteronuclear control demonstration inspired by the low-power, offset-tolerant (^{14}mathrm N)-decoupling work of Nehra, Agarwal, and Nishiyama. It represents one observed (^{1}mathrm H) spin coupled to one controlled (^{14}mathrm N) spin, with an effective scalar coupling of 1500 and an 18.8 T field (approximately 800 MHz (^{1}mathrm H)). The model uses an effective interaction rather than a full quadrupolar/MAS Hamiltonian; its source describes the quadrupolar/MAS response as compressed into an effective nitrogen-offset ensemble.

Only the (^{14}mathrm N) channel is controlled, while normalized transverse (^{1}mathrm H) states are the preservation targets. A phase-only waveform is optimized with GRAPE/`fmaxnewton` from an XY-type phase seed, with Bloch–Siegert correction enabled on the driven nitrogen channel. The 120 elements are 10 μs each (1.2 ms total), at a nominal 20 kHz RF field. The comparison is a constant-phase CW-like waveform with the same duration and power—not a separate optimized control. Training uses seven offsets from −12 to +12 kHz and three (B_1) scales (0.95, 1.00, 1.05); evaluation uses 61 offsets from −20 to +20 kHz and nine scales from 0.90 to 1.10.

The source plots the phase cycle, offset- and (B_1)-averaged proton-preservation fidelities, and offset/(B_1) maps; it also prints training-ensemble mean and worst-case fidelities for the two waveforms. Those are outputs the script is designed to calculate, not results asserted by this note. The source names Nehra, Agarwal, and Nishiyama but gives no complete bibliographic reference; no missing citation details are invented.
