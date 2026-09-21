# examples/nmr_proteins/hcanh_gb1.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_proteins/hcanh_gb1.m`
- Signature: `hcanh_gb1()`
- Total lines: 79

## Purpose

Simulated H(CA)NH spectrum of GB1 protein. It is assumed that only the backbone is 13C,15N-labelled. Calculation time: minutes, faster with a Tesla A100 GPU.

## Physical / mathematical content

- Protein NMR examples. These files specialise liquid-state pulse sequences to labelled biomolecules, exploiting one-bond and two-bond heteronuclear couplings, coherence pathway filtering, selective decoupling, and high-dimensional indirect detection.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Implementation structure

- Simulated H(CA)NH spectrum of GB1 protein. It is assumed that
- only the backbone is 13C,15N-labelled.
- Calculation time: minutes, faster with a Tesla A100 GPU.
- Protein data import
- Magnet field
- Tolerances
- Basis set
- Algorithmic options
- Spinach housekeeping
- Sequence parameters
- Simulation
- Apodisation

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `protein()`, `create()`, `basis()`, `liquid()`, `apodisation()`, `fftshift()`, `conj()`, `kfigure()`, `plot_3d()`.
