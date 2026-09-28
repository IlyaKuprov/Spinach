# tests/kernel/test_dynamic_examples_smoke.m

- Signature: `result=test_dynamic_examples_smoke()`

## Purpose

Smoke-tests compact liquid-state NMR example paths through acquisition, spectrum processing, and plotting. Returns a regression test result with explanatory messages.

## Physical / mathematical content

- The one-spin `1H` system has a 14.1 T field and zero Zeeman offset. Its transverse free-induction decay (FID) is checked against a constant value of `0.5` at all eight points; its eight-point Fourier transform is checked against a spectrum with only the centred DC point nonzero.
- The two-spin `1H` CT-COSY system has a 5.9 T field, scalar Zeeman values of `1.00` and `3.00`, and a scalar coupling of `7.0` between the spins.

## Numerical / algorithmic content

- Both systems use the `sphten-liouv` formalism with `approximation='none'` and run through `liquid(...,'nmr')`.
- The one-dimensional path calls `@acquire` with offset `0`, sweep `1000`, eight acquired points, eight spectrum points, and Hz axis units. It computes `fftshift(fft(fid,parameters.zerofill))` and compares the FID and spectrum with their references using absolute and relative tolerances of `1e-12`.
- The CT-COSY path calls `@ct_cosy` with offset `500`, sweep `[2000 2000]`, an `[8 8]` FID, a `[16 16]` spectrum, and Hz axis units. It computes `fftshift(fft2(fid,parameters.zerofill(2),parameters.zerofill(1)))` and checks the FID and spectrum sizes. It also checks the FID norm (`5.6305493030431162e+00`), spectrum norm (`9.0088788848689859e+01`), absolute-spectrum sum (`6.3367922545678312e+02`), and maximum (`3.1850510896725737e+01`) against deterministic references.

## Outputs

- `result` — regression test result with explanatory messages. The test also checks that `plot_1d` produces line Y-data equal to the real one-spin spectrum, and that `plot_2d` returns the transposed absolute CT-COSY spectrum, axes matching the spectrum dimensions, and one contour object.

## Implementation structure

- Announces the dynamic plotting smoke test and creates a result under `examples/dynamic_examples_smoke`.
- Saves the default figure visibility, makes figures invisible, and runs the one-dimensional acquisition and two-dimensional CT-COSY checks in separate local functions. Each plotting check closes its figure.
- An `onCleanup` handler closes all figures and restores the saved default figure visibility after success or failure.

Attribution: ilya.kuprov@weizmann.ac.il