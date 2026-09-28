# examples/nmr_spen/ufcosy_2spin.m

- Signature: `ufcosy_2spin()`

## Purpose

Ultrafast COSY for a coupled two-spin system. Calculation time: minutes on NVidia Tesla A100, much longer on CPU Jean-Nicolas Dumez Ludmilla Guduff

## Physical / mathematical content

- The source defines two 1H spins at 14.0 T with a 10 Hz scalar coupling. It models a 15 mm sample with 500 spatial points; the diffusion coefficient and flow are set to zero.
- The imaging simulation uses `@spencosy` with encoding and coherence selection. It acquires 512 points over 128 loops at 0.5 μs dwell; sequence settings include `pulsenpoints=1000`, `nWURST=40`, `Te=15 ms`, and `BW=10 kHz`, with `Ga=0.50`, `Ge=0.01`, and `Gp=0.47 T/m` and `Tp=1 ms`.

## Numerical / algorithmic content

- The script Fourier-transforms the FID along its second dimension, then displays the magnitude as a contour plot.

## Implementation structure

- Builds the spin system and basis, configures the sample and sequence, runs imaging, and processes and plots the result.
