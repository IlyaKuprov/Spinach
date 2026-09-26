# examples/nmr_proteins/methyl_noesy_gb1.m

- Signature: `methyl_noesy_gb1()`

## Purpose

1H-1H NOESY spectrum of GB1 with everything deuterated except methyl groups. Deuteria are kept in the spin system because they are a part of the coupling network; methyl group rotati- on is not accounted for in this simulation. Calculation time: hours.

## Physical / mathematical content

- Protein NMR examples. These files specialise liquid-state pulse sequences to labelled biomolecules, exploiting one-bond and two-bond heteronuclear couplings, coherence pathway filtering, selective decoupling, and high-dimensional indirect detection.
- The relaxation model is Redfield-type perturbation theory: fluctuating interactions enter through correlation functions or spectral densities and generate a linear relaxation superoperator.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Implementation structure

- 1H-1H NOESY spectrum of GB1 with everything deuterated except
- methyl groups. Deuteria are kept in the spin system because
- they are a part of the coupling network; methyl group rotati-
- on is not accounted for in this simulation.
- Calculation time: hours.
- Protein data import
- Magnet field
- Tolerances
- Relaxation theory
- Basis set
- Algorithmic options
- Create the spin system structure
