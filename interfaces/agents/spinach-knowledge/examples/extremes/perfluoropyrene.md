# examples/extremes/perfluoropyrene.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/extremes/perfluoropyrene.m`
- Signature: `perfluoropyrene()`
- Total lines: 64

## Purpose

X-band pulsed ESR spectrum of perfluoropyrene cation radical, computed using brute force operator algebra in the full 4,194,304 -dimensional Liouville space. This is deliberate -a much faster calculation is, of course, possible with a restricted basis set. This calculation requires at least 64GB of RAM and illustrates the per- formance of trajectory-level state space restriction in Spinach. Calculation time: minutes

## Physical / mathematical content

- Extreme-regime examples. These scripts exercise Spinach in unusually large, stiff, high-field, low-field, or otherwise numerically demanding regimes where approximations, conditioning, and basis-size control are central.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

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

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `g2spinach()`, `gparse()`, `create()`, `basis()`, `state()`, `liquid()`, `apodisation()`, `fftshift()`, `kfigure()`, `plot_1d()`.
