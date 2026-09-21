# examples/kinetics/mas_exchange_1.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/kinetics/mas_exchange_1.m`
- Signature: `mas_exchange_1()`
- Total lines: 69

## Purpose

Two-site position exchange for a deuterium nucleus. The sites differ in the chemical shift and the orientation of the quad- rupolar tensor. Calculation time: seconds.

## Physical / mathematical content

- Chemical-kinetics examples. The files couple spin dynamics to exchange, pumping, or nonlinear reaction networks represented by kinetic generators in Liouville space.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.
- Quadrupolar physics is relevant: nuclei with spin > 1/2 interact with the electric field gradient tensor, introducing second-rank anisotropy, asymmetry, and overtone or MQ phenomena.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Implementation structure

- Two-site position exchange for a deuterium nucleus. The sites
- differ in the chemical shift and the orientation of the quad-
- rupolar tensor.
- Calculation time: seconds.
- System specification
- Spin system
- Quadrupolar interactions
- Chemical shifts
- Chemical exchange
- Basis set
- Spinach housekeeping
- Sequence parameters

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `weblab2nqi()`, `acos()`, `create()`, `basis()`, `state()`, `singlerot()`, `apodisation()`, `fftshift()`, `kfigure()`, `plot_1d()`.
