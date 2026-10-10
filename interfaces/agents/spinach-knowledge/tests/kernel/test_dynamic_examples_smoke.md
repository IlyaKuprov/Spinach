# tests/kernel/test_dynamic_examples_smoke.m

**Source:** [https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_dynamic_examples_smoke.m](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_dynamic_examples_smoke.m)

## Purpose

Regression test that exercises compact dynamic example-stage calculations with plotting. It runs short liquid-state NMR calculations adapted from plotting examples, processes the deterministic signals, and verifies the plotted graphics objects under invisible offscreen figures.

## What the test establishes

- **Analytic one-spin acquisition:** a zero-offset `1H` example yields an eight-point FID constant at `0.5`. Its eight-point Fourier spectrum has only the centred DC bin, of amplitude `8×0.5`. Both data arrays and the one-dimensional plot output are compared with the analytic references at `1e-12` absolute and relative tolerance.
- **Two-spin CT-COSY regression:** a compact pair of `1H` spins at the original 5.9 T proton frequency with 7 Hz coupling, 500 Hz offset and `[2000 2000]` Hz sweeps produces an `8×8` FID and a `16×16` zero-filled spectrum. The reference FID and spectrum 2-norms are `5.6305493030431162e+00` and `9.0088788848689859e+01` (`1e-10` tolerances); the absolute-spectrum sum is `6.3367922545678312e+02` (`1e-9` absolute, `1e-10` relative) and its maximum `3.1850510896725737e+01` (`1e-10`). These are fixed regression targets, not newly measured results.
- **Offscreen rendering:** invisible figure creation must not alter the underlying numerical results. The two-dimensional plot returns the transpose of the absolute spectrum (`1e-12`), both frequency axes have the appropriate zero-filled lengths, and one contour object is present.

The one-spin and CT-COSY proton carriers are fixed at `3772062842.904` and `1578380905.8960001` rad/s, respectively, by dividing those angular frequencies by `spin('1H')` to obtain the field. Updating magnetic moments therefore does not change the regression's physical frequencies or require relaxed numerical tolerances.

## Inputs and outputs

```matlab
result=test_dynamic_examples_smoke()
```

- **Output:** `result` — regression test result with explanatory messages, accumulated through `test_close` and `test_true` assertions.
- Takes no inputs.

## References

- [Spinach — tests/kernel/test_dynamic_examples_smoke.m (GitHub)](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_dynamic_examples_smoke.m)
