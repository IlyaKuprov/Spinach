# examples/parahydrogen/ortho_deuterium.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/parahydrogen/ortho_deuterium.m`
- Signature: `ortho_deuterium()`
- Total lines: 68

## Purpose

Ortho-deuteration simulation for acrylonitrile in Figure 1 of the paper by Natterer, Greve, and Bargon: Simulation time: seconds

## Physical / mathematical content

- Parahydrogen examples. The physical motif is highly non-Boltzmann singlet order imported from para-H2 and converted into observable nuclear magnetisation through hydrogenation, exchange, or catalytic transfer processes.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.
- The relevant state manifold is the singlet/triplet decomposition, where permutation symmetry controls selection rules, relaxation susceptibility, and convertibility to ordinary magnetisation.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Implementation structure

- Ortho-deuteration simulation for acrylonitrile in Figure 1 of
- the paper by Natterer, Greve, and Bargon:
- Simulation time: seconds
- Bargon's magnet
- Deuteration product
- Hilbert space
- Spinach housekeeping
- Continuous deuteration
- Singlet and quintet on deuterium
- Experiment parameters
- Simulation
- Apodisation and sign flip

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `deut_pair()`, `state()`, `operator()`, `liquid()`, `apodisation()`, `fftshift()`, `kfigure()`, `plot_1d()`.
