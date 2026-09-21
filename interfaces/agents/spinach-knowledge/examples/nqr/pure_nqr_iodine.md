# examples/nqr/pure_nqr_iodine.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nqr/pure_nqr_iodine.m`
- Signature: `pure_nqr_iodine()`
- Total lines: 53

## Purpose

Powder NQR spectrum of a system with a single 127I nucleus. Calculation time: seconds

## Physical / mathematical content

- NQR examples. The Hamiltonian is dominated by quadrupolar interaction with little or no Zeeman field, so transition frequencies reflect electric field gradients and asymmetry parameters.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Implementation structure

- Powder NQR spectrum of a system with a single 127I nucleus.
- Calculation time: seconds
- System specification
- Formalism and basis
- Relaxation theory
- Spinach housekeeping
- Experiment parameters
- Simulation
- Apodisation
- Fourier transform
- Plotting

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `eeqq2nqi()`, `create()`, `basis()`, `state()`, `operator()`, `powder()`, `apodisation()`, `fftshift()`, `kfigure()`, `plot_1d()`.
