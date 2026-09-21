# examples/nmr_proteins/hncoca_simple.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_proteins/hncoca_simple.m`
- Signature: `hncoca_simple()`
- Total lines: 76

## Purpose

A minimal HNCOCA pulse sequence simulation. Calculation time: seconds.

## Physical / mathematical content

- Protein NMR examples. These files specialise liquid-state pulse sequences to labelled biomolecules, exploiting one-bond and two-bond heteronuclear couplings, coherence pathway filtering, selective decoupling, and high-dimensional indirect detection.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Implementation structure

- A minimal HNCOCA pulse sequence simulation.
- Calculation time: seconds.
- Magnet field
- Spin system
- Interactions
- Basis set
- Spinach housekeeping
- Sequence parameters
- Simulation
- Apodisation
- F3 Fourier transform
- Absorption part of F3 signal

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `liquid()`, `apodisation()`, `fftshift()`, `conj()`, `kfigure()`, `plot_3d()`.
