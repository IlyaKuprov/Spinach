# examples/nmr_proteins/expt_data/hsqc_ubiquitin_expt.m

- Signature: `hsqc_ubiquitin_expt()`

## Purpose

Experimental HSQC spectrum of human ubiquitin. Donghan Lee (Max Planck Institute) Ilya Kuprov (University of Southampton)

## Physical / mathematical content

- Protein NMR examples. These files specialise liquid-state pulse sequences to labelled biomolecules, exploiting one-bond and two-bond heteronuclear couplings, coherence pathway filtering, selective decoupling, and high-dimensional indirect detection.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Implementation structure

- Experimental HSQC spectrum of human ubiquitin.
- Donghan Lee (Max Planck Institute)
- Ilya Kuprov (University of Southampton)
- Load the file
- Phasing
- Apodisation
- F2 Fourier transform
- Form States signal
- F1 Fourier transform
- Plotting
