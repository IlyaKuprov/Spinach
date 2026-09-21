# examples/nmr_proteins/hcch_cosy_simple.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_proteins/hcch_cosy_simple.m`
- Signature: `hcch_cosy_simple()`
- Total lines: 81

## Purpose

3D HCCH COSY experiment on a small protein fragment. Calculation time: seconds.

## Physical / mathematical content

- Protein NMR examples. These files specialise liquid-state pulse sequences to labelled biomolecules, exploiting one-bond and two-bond heteronuclear couplings, coherence pathway filtering, selective decoupling, and high-dimensional indirect detection.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Implementation structure

- 3D HCCH COSY experiment on a small protein fragment.
- Calculation time: seconds.
- Spin system
- Interactions
- Basis set
- Sequence parameters
- Spinach housekeeping
- Simulation
- Apodisation
- F3 Fourier transform
- Absorption part of F3 signal
- F2 Fourier transform

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `liquid()`, `apodisation()`, `fftshift()`, `conj()`, `kfigure()`, `plot_3d()`.
