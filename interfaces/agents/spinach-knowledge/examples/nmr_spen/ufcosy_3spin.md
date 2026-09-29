# examples/nmr_spen/ufcosy_3spin.m

- Function: `ufcosy_3spin()`
- Source: [examples/nmr_spen/ufcosy_3spin.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_spen/ufcosy_3spin.m)

This is an ultrafast COSY simulation for three coupled 1H spins, passed to Spinach imaging as the `spencosy` sequence. The field is 14.0 T. The chemical-shift inputs are `{3.70, 3.92, 4.50} - {4.0, 4.0, 4.0}` ppm, or -0.30, -0.08, and +0.50 ppm relative to 4.0 ppm. Pairwise scalar couplings are 10 Hz (spins 1-2), 12 Hz (2-3), and 4 Hz (1-3). The file is a parameterised spin simulation, not an imported experimental spectrum.

The sample model is 15 mm long with 500 spatial points and period-3 spatial differentiation. Flow and diffusion are zero. Uniform initial and receive phantoms use the 1H longitudinal state and transverse receive state. Relaxation phantoms and relaxation operators are both empty in this example; it therefore does not add a relaxation model through these parameters.

The sequence inputs match the two-spin demonstration: 1H acquisition at zero offset, 0.5 microsecond dwell, 512 points and 128 loops, with acquisition gradient `Ga=0.50 T/m`. Encoding uses 1000 pulse points, `nWURST=40`, `Te=15 ms`, 10 kHz bandwidth, and `Ge=0.01 T/m`; coherence selection uses `Gp=0.47 T/m` for `Tp=1 ms`. The wrapper provides these encoding and selection parameters to `spencosy`, but does not itself specify a full pulse/gradient waveform or detailed coherence-pathway schedule.

The simulated FID is Fourier-transformed along its second dimension, shifted, and shown as a magnitude contour plot. The basis is sphten-liouv with no approximation; PT is disabled and the greedy algorithm enabled. There is no measured-data comparison or accuracy bound in the source. Its machine-time comment is not reported as a validated runtime.
