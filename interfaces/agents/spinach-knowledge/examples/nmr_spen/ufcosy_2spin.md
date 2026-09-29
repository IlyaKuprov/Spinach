# examples/nmr_spen/ufcosy_2spin.m

- Function: `ufcosy_2spin()`
- Source: [examples/nmr_spen/ufcosy_2spin.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_spen/ufcosy_2spin.m)

This example simulates ultrafast COSY for a coupled two-proton system, using the `spencosy` imaging sequence. Both spins are 1H at 14.0 T. The chemical-shift inputs are `{3.70, 4.50} - {4.0, 4.0}` ppm, i.e. offsets of -0.30 and +0.50 ppm from the shared 4.0 ppm reference. Their scalar coupling is 10 Hz. This is a specified spin model, not a dataset imported from a measurement.

The spatial model is a 15 mm sample on 500 points, with the period-3 derivative setting. Flow and diffusion are set to zero. Initial and receive phantoms are uniform, with 1H longitudinal initial state and 1H transverse detection state. Relaxation uses a relaxation operator and zero-valued relaxation phantom.

The acquisition is 1H-selective, with zero offset, 0.5 microsecond dwell, 512 points per loop and 128 loops; the acquisition gradient `Ga` is 0.50 T/m. Encoding inputs are 1000 pulse points, `nWURST=40`, `Te=15 ms`, bandwidth 10 kHz, and encoding gradient `Ge=0.01 T/m`. Coherence selection uses `Gp=0.47 T/m` for `Tp=1 ms`. These are the values supplied to `spencosy`; the example does not define a complete time-resolved waveform schedule in this wrapper, so the internal chirp and gradient timing is not elaborated beyond those inputs.

After imaging, the code Fourier-transforms the FID along dimension 2, shifts the spectrum, and plots its magnitude as contours. The basis uses the sphten-liouv formalism without an approximation; the example disables the PT option and enables the greedy algorithm. No measured spectrum, fit, or accuracy limit is supplied, and the source's A100/CPU timing comment is not a validated runtime result.
