# examples/nmr_proteins/expt_data/hnco_ubiquitin_expt.m

- Signature: `hnco_ubiquitin_expt()`

## Purpose

Experimental HNCO spectrum of human ubiquitin. Donghan Lee (Max Planck Institute) Ilya Kuprov (University of Southampton)

## Physical / mathematical content

- Protein NMR examples. These files specialise liquid-state pulse sequences to labelled biomolecules, exploiting one-bond and two-bond heteronuclear couplings, coherence pathway filtering, selective decoupling, and high-dimensional indirect detection.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Implementation structure

- Experimental HNCO spectrum of human ubiquitin.
- Donghan Lee (Max Planck Institute)
- Ilya Kuprov (University of Southampton)
- Data loading and truncation
- Apodisation
- F3 processing
- F2 processing
- F1 processing
- Window shifting
- Baseline correction
- Water signal elimination
- Parameters
