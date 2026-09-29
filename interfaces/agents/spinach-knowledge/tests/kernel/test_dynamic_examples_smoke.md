# tests/kernel/test_dynamic_examples_smoke.m

**Source:** [https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_dynamic_examples_smoke.m](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_dynamic_examples_smoke.m)

## Purpose

Regression test that exercises compact dynamic example-stage calculations with plotting. It runs short liquid-state NMR calculations adapted from plotting examples, processes the deterministic signals, and verifies the plotted graphics objects under invisible offscreen figures.

## Behavior

- Announces the test target with `TESTING: Dynamic plotting example smoke paths` and initializes a test result named `examples/dynamic_examples_smoke` with the description `Dynamic plotting example smoke paths` and the requirement that compact example-stage calculations must run, process, and plot deterministic spectra.
- Forces invisible figures by saving and overriding `groot`'s `defaultFigureVisible` to `'off'`; an `onCleanup` object restores the original visibility and closes all figures with `close all force` after the test, whether it succeeds or fails.
- Runs two subtests:
  - **One-dimensional acquisition path** (`local_test_acquire_1d`):
    - Builds a zero-offset one-spin Liouville-space system with `sys.magnet=14.1`, `sys.isotopes={'1H'}`, `inter.zeeman.scalar={0}`, `bas.formalism='sphten-liouv'`, `bas.approximation='none'` via `test_spin_system`.
    - Sets up a compact free-induction acquisition: `parameters.spins={'1H'}`, `rho0` and `coil` both `state(spin_system,'L+','1H')`, `decouple={}`, `offset=0`, `sweep=1000`, `npoints=8`, `zerofill=8`, `axis_units='Hz'`, `invert_axis=0`.
    - Runs `liquid(spin_system,@acquire,parameters,'nmr')` and computes `spectrum=fftshift(fft(fid,parameters.zerofill))`.
    - Verifies against analytic references: the zero-offset FID is constant at `0.5` for all 8 points, and its Fourier transform has only the centred DC point, with bin `zerofill/2+1` equal to `zerofill*fid_ref(1)`; both comparisons use tolerances `1e-12`.
    - Plots the real part of the spectrum with `plot_1d(spin_system,real(spectrum),parameters,'k-')` on an invisible `kfigure`, retrieves the first line object's `YData`, and checks it matches `real(spectrum)` to `1e-12`.
  - **Two-dimensional CT-COSY path** (`local_test_ct_cosy_2d`):
    - Builds a two-spin system with `sys.isotopes={'1H','1H'}`, `sys.magnet=5.9`, `inter.zeeman.scalar={1.00 3.00}`, `inter.coupling.scalar{1,2}=7.0`, `inter.coupling.scalar{2,2}=0`, `bas.formalism='sphten-liouv'`, `bas.approximation='none'`.
    - Uses compact point counts: `offset=500`, `sweep=[2000 2000]`, `npoints=[8 8]`, `zerofill=[16 16]`, `spins={'1H'}`, `axis_units='Hz'`.
    - Runs `liquid(spin_system,@ct_cosy,parameters,'nmr')` and computes `spectrum=fftshift(fft2(fid,parameters.zerofill(2),parameters.zerofill(1)))`, then `abs_spectrum=abs(spectrum)`.
    - Checks deterministic sizes and invariants:
      - FID size equals `parameters.npoints` and spectrum size equals `parameters.zerofill`.
      - FID 2-norm equals `5.6305493030431162e+00` (tolerances `1e-10`).
      - Spectrum 2-norm equals `9.0088788848689859e+01` (tolerances `1e-10`).
      - Absolute spectrum sum equals `6.3367922545678312e+02` (tolerances `1e-9` and `1e-10`).
      - Absolute spectrum maximum equals `3.1850510896725737e+01` (tolerances `1e-10`).
    - Plots the processed spectrum with `scale_figure([1.5 2.0])` and `plot_2d(spin_system,abs_spectrum,parameters,6,[0.10 0.50 0.10 0.50],2,64,4,'positive')` on an invisible `kfigure`.
    - Verifies the returned `plot_spectrum` matches `transpose(abs_spectrum)` to `1e-12`, that `axis_f1` and `axis_f2` lengths equal `zerofill(2)` and `zerofill(1)` respectively, and that exactly one contour object exists on the figure.
- Closes each figure after its checks.

## Inputs and outputs

```matlab
result=test_dynamic_examples_smoke()
```

- **Output:** `result` — regression test result with explanatory messages, accumulated through `test_close` and `test_true` assertions.
- Takes no inputs.

## References

- [Spinach — tests/kernel/test_dynamic_examples_smoke.m (GitHub)](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_dynamic_examples_smoke.m)
