# examples/nmr_spen/ufmq_6spin.m

- Signature: `ufmq_6spin()`
- Source: [`examples/nmr_spen/ufmq_6spin.m`](https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_spen/ufmq_6spin.m)

## Experiment and spin model

This example simulates a six-quantum (6Q), ultrafast MaxQ NMR experiment on six coupled protons, with spatial encoding and diffusion. It is a Spinach simulation, not an imported experimental spectrum. The source credits Maria Grazia Concilio, Ilya Kuprov and Jean-Nicolas Dumez. The 14.1 T model has six 1H sites with shifts −1.0, −0.5, 0, +0.3, +0.7, and +0.9 ppm. The listed 3J couplings are 8.00 Hz for pairs 1–2, 2–3, 3–4, and 4–5; the 4J entries are 4.00 Hz for 1–3, 1–4, 2–4, 3–5, 1–6, and 5–6; the 5J entries are 4.00 Hz for 2–5, 2–6, and 4–6; and the 6J entries are 2.00 Hz for 1–5 and 3–6. The source sets coherence order to +6 and uses the full sphten-liouv basis without approximation.

## Spatial encoding and acquisition

The imaging phantom is 0.015 m long with 300 points, uniform initial and receive phantoms, zero flow, and diffusion coefficient 18 × 10⁻¹⁰ m²/s. The initial state is 1H Lz and the receiver is 1H L+. Acquisition uses 1H, zero offset, 6 μs dwell, 120 points, and 50 loops. The maximum k value is calculated as npoints/dims; the acquisition-gradient amplitude is calculated from that k value and the 720 μs acquisition duration. The WURST encoding parameters are 500 pulse points, 40 WURST cycles, Te = 15 ms, bandwidth = 15 kHz, and encoding gradient Ge = 0.023 T/m, with a 41 ms delay. No field-drop schedule is specified in this source.

## Output and limits

The script calls `imaging(...,@ufmq,...)`, displays the real k-space echo array, Fourier transforms along the conventional dimension, and plots the magnitude spectrum in ppm. The source comments estimate hours on an NVIDIA Tesla A100 and much longer on CPU; this is a source estimate, not a measured runtime. The GPU-enable setting is commented out; the active algorithm options disable `pt` and enable `zte` and `greedy`.
