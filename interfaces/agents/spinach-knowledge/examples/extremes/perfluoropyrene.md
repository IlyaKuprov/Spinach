# examples/extremes/perfluoropyrene.m

- Signature: `perfluoropyrene()`

## Purpose

X-band pulsed ESR spectrum of perfluoropyrene cation radical, computed using brute-force operator algebra in the full 4,194,304-dimensional Liouville space. This is deliberate—a much faster calculation is possible with a restricted basis set. The calculation requires at least 64 GB of RAM and illustrates the performance of trajectory-level state-space restriction in Spinach. Calculation time: minutes

## Physical / mathematical content

- The example computes the X-band pulsed ESR spectrum of the perfluoropyrene cation radical in a 4,194,304-dimensional Liouville space; the full-space calculation is retained deliberately to demonstrate trajectory-level state-space restriction.
- The detected time-domain signal is apodised and Fourier-transformed to produce the ESR spectrum.

## Numerical / algorithmic content

- The script propagates the pulsed ESR signal in the full Liouville space, applies apodisation, and uses an FFT for the plotted spectrum.

## Implementation structure

- X-band pulsed ESR spectrum of perfluoropyrene cation radical, computed
- using brute force operator algebra in the full 4,194,304 -dimensional
- Liouville space. This is deliberate -a much faster calculation is, of
- course, possible with a restricted basis set.
- This calculation requires at least 64GB of RAM and illustrates the per-
- formance of trajectory-level state space restriction in Spinach.
- Calculation time: minutes
- Ignore coordinate information (HFCs provided)
- Read the spin system properties (vacuum DFT calculation)
- Magnet induction
- Relaxation theory
- Basis set
