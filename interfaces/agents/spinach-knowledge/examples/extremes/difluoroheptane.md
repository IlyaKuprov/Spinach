# examples/extremes/difluoroheptane.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/extremes/difluoroheptane.m`
- Signature: `difluoroheptane()`
- Total lines: 122

## Purpose

19F NMR spectrum of anti-3,4-difluoroheptane (16 spins) by explicit time-domain evolution in Liouville space. WARNING: needs 32 CPU cores, 128 GB of RAM and a Titan V or later. Run time on the above: minutes

## Physical / mathematical content

- Extreme-regime examples. These scripts exercise Spinach in unusually large, stiff, high-field, low-field, or otherwise numerically demanding regimes where approximations, conditioning, and basis-size control are central.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Implementation structure

- 19F NMR spectrum of anti-3,4-difluoroheptane (16 spins) by
- explicit time-domain evolution in Liouville space.
- WARNING: needs 32 CPU cores, 128 GB of RAM and
- a Titan V or later.
- Run time on the above: minutes
- Magnet induction
- Isotopes
- Shifts
- J-couplings
- Basis set
- Greedy parallelisation
- Spinach housekeeping

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `false()`, `create()`, `basis()`, `state()`, `liquid()`, `apodisation()`, `fftshift()`, `kfigure()`, `plot_1d()`.
