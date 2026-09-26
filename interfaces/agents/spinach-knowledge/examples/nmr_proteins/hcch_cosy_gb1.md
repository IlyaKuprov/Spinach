# examples/nmr_proteins/hcch_cosy_gb1.m

- Signature: `hcch_cosy_gb1()`

## Purpose

3D HCCH COSY experiment on GB1 protein. Calculation time: hours.

## Physical / mathematical content

- Protein NMR examples. These files specialise liquid-state pulse sequences to labelled biomolecules, exploiting one-bond and two-bond heteronuclear couplings, coherence pathway filtering, selective decoupling, and high-dimensional indirect detection.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Implementation structure

- 3D HCCH COSY experiment on GB1 protein.
- Calculation time: hours.
- Protein data import
- Magnet field
- Tolerances
- Basis set
- Algorithmic options
- Sequence parameters
- Create the spin system structure
- Kill nitrogens (not relevant)
- Build the basis
- Simulation
