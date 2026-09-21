# examples/nmr_liquids/inv_rec_strychnine.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_liquids/inv_rec_strychnine.m`
- Signature: `inv_rec_strychnine()`
- Total lines: 62

## Purpose

1H inversion-recovery experiment on strychnine at 250 MHz. Calculation time: minutes

## Physical / mathematical content

- Liquid-state NMR examples. The physics is scalar-coupling-mediated coherence transfer in weakly or moderately coupled spin systems, often in Liouville space. Typical mechanisms include INEPT-style polarisation transfer, J-refocusing, phase cycling, indirect evolution, and multidimensional detection.
- The relaxation model is Redfield-type perturbation theory: fluctuating interactions enter through correlation functions or spectral densities and generate a linear relaxation superoperator.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Implementation structure

- 1H inversion-recovery experiment on strychnine at 250 MHz.
- Calculation time: minutes
- Read spin system properties
- Magnetic induction
- Maximum distance to consider
- Greedy parallelisation
- Relaxation theory parameters
- Basis set
- Spinach housekeeping
- Sequence parameters
- Simulation
- Apodisation

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `strychnine()`, `create()`, `basis()`, `liquid()`, `apodisation()`, `fftshift()`, `kfigure()`, `plot_1d()`.
