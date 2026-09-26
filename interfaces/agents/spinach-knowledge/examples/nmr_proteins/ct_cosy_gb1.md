# examples/nmr_proteins/ct_cosy_gb1.m

- Signature: `ct_cosy_gb1()`

## Purpose

Constant-time COSY experiment simulation for the GB1 protein. Simulation time: minutes, faster with a Tesla A100 GPU.

## Physical / mathematical content

- Protein NMR examples. These files specialise liquid-state pulse sequences to labelled biomolecules, exploiting one-bond and two-bond heteronuclear couplings, coherence pathway filtering, selective decoupling, and high-dimensional indirect detection.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Implementation structure

- Constant-time COSY experiment simulation for the
- GB1 protein.
- Simulation time: minutes, faster with a Tesla A100 GPU.
- Protein data import
- Magnet field
- Tolerances
- Basis set
- Algorithmic options
- Sequence parameters
- Spinach housekeeping
- Kill carbons and nitrogens (protein assumed unlabelled)
- Build the basis
