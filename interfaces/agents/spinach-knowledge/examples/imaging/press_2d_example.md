# examples/imaging/press_2d_example.m

- Signature: `press_2d_example()`

## Purpose

Demonstrates 2D PRESS localisation and spectral readout for three spin-pair components placed at separate positions in a 2D sample. Runtime is estimated in minutes, faster with a Tesla V100 GPU.

## Spin systems and sample

The 3.0 T model has six protons in three scalar-coupled pairs, with chemical-shift values `{-3,+3}`, `{-2,+2}`, and `{-1,+1}`, and pair couplings 10, 20, and 30. It uses the `sphten-liouv` formalism, `IK-2` approximation, proximity level 1, and scalar-coupling connectivity. The sample dimensions are `[0.30 0.25] m` on a `[108 90]` grid. Three circular density masks are centred at grid coordinates (15,15), (30,40), and (70,80), each with radius 10 points.

## PRESS acquisition and output

The image size is `[101 105]`; the selection-gradient amplitudes are `[25 25] mT/m`. Two RF pulses have phase `pi/2`, amplitude `2*pi*5000`, and durations `0.5e-4` and `1.0e-4 s`; the executed frequency list is `{-120e3 -100e3}` and maximum ranks are `{3 3}`. The source comments list scan settings of {-120,-100} kHz for A, {-80,-10} kHz for B, and {+30,+100} kHz for C.

The example plots the combined sample phantom and active-volume diagnostic, runs the 2D PRESS sequence, applies square-cosine apodisation, and plots the magnitude Fourier-transformed spectrum.
