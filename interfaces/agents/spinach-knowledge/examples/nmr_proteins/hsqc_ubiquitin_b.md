# examples/nmr_proteins/hsqc_ubiquitin_b.m

- Signature: `hsqc_ubiquitin_b()`

## Purpose

1H-15N HSQC of human ubiquitin, without 1H decoupling in F1 and 15N decoupling in F2. Nitrogen-proton multiplicity is retained in both dimensions. Calculation time: hours, faster with a Tesla A100 GPU. Zenawi Welderufael Luke Edwards Ilya Kuprov

## Physical / mathematical content

- Protein NMR examples. These files specialise liquid-state pulse sequences to labelled biomolecules, exploiting one-bond and two-bond heteronuclear couplings, coherence pathway filtering, selective decoupling, and high-dimensional indirect detection.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Implementation structure

- 1H-15N HSQC of human ubiquitin, without 1H decoupling in F1 and
- 15N decoupling in F2. Nitrogen-proton multiplicity is retained
- in both dimensions.
- Calculation time: hours, faster with a Tesla A100 GPU.
- Zenawi Welderufael
- Luke Edwards
- Ilya Kuprov
- Protein data import
- Magnet field
- Tolerances
- Basis set
- Algorithmic options
