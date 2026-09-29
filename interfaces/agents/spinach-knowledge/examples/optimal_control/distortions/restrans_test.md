# examples/optimal_control/distortions/restrans_test.m

- MATLAB implementation: [examples/optimal_control/distortions/restrans_test.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/optimal_control/distortions/restrans_test.m)

Source: [examples/optimal_control/distortions/restrans_test.m](../../../../../../examples/optimal_control/distortions/restrans_test.m)

- Signature: `restrans_test()`

## Purpose

Demonstrates the time-domain resonator transform on a generated two-channel waveform for two resonance settings. The pulse is random Gaussian data, not a square pulse, and the script does not formulate an optimal-control problem.

## Input waveform and model settings

The script draws a 51-by-2 array from MATLAB `randn`, then zeros the first and last three samples in both columns. The columns are passed as X and Y pulse components. Each sample is assigned a 1 μs slice duration, and both calls use Q = 50 with the piecewise-constant model option `'pwc'` and 100 model points.

The two side-by-side calls use resonance angular frequencies of 2π × 600 MHz and 2π × 60 MHz, respectively, and are titled as ¹H and ¹⁵N at 14.1 T. The source labels the example as an RLC response simulation. No measured input waveform or experimental response file is loaded.

## Output and limits

The script delegates each calculation and its time-domain response display to `restrans`, placing the two cases in separate subplots. The random waveform is not seeded in the script, so its values can differ between runs. This is an illustrative model calculation; the example does not report a hardware comparison or a quantitative fit to measured data.
