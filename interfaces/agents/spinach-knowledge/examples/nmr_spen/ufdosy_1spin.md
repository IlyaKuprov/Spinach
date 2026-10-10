# examples/nmr_spen/ufdosy_1spin.m

- Function: `ufdosy_1spin()`
- Source: [examples/nmr_spen/ufdosy_1spin.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_spen/ufdosy_1spin.m)

This example simulates ultrafast DOSY for one 1H spin using the `spendosy` sequence. The field is 14.1 T and the shift input is 7.0 ppm. It is a one-spin simulated phantom, not an imported measurement or a fitted diffusion result.

The sample is 15 mm long, discretised at 3000 points with the period-7 spatial derivative setting. Flow is zero and the diffusion coefficient is `8e-10 m^2/s`. The source invokes the NMR assumption, uses uniform initial and receive phantoms, and sets longitudinal initial magnetisation and transverse 1H detection. Relaxation phantoms and relaxation operators are empty.

Acquisition inputs are 1H, `td=100 ms`, dwell 1.5 microseconds, 128 points, 256 loops, acquisition gradient `Ga=0.52 T/m`, and offset 3600 Hz. The source does not explain `td` further, so its physical interpretation is not expanded here. Encoding uses 1000 pulse points, smoothing factor 0.1, `Te=1.5 ms`, `Tau=1.6 ms`, bandwidth 110 kHz, encoding gradient `Ge=0.2535 T/m`, and the smoothed chirp type. These are sequence inputs; the wrapper does not state a more detailed gradient/chirp waveform schedule.

The code runs imaging with `@spendosy`, Fourier-transforms both FID dimensions, and plots the magnitude. It computes the chemical-shift axis from the offset, acquisition width, and magnet frequency, and the displacement axis from the dwell and acquisition gradient; the latter is displayed in millimetres. This gives a shift-versus-position image for the specified diffusion simulation. No experimental comparison, diffusion-fit validation, or accuracy limit is included. The source's A100/CPU time comment is not a validated runtime.
