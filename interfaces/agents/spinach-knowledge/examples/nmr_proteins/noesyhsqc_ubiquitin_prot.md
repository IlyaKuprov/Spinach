# examples/nmr_proteins/noesyhsqc_ubiquitin_prot.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_proteins/noesyhsqc_ubiquitin_prot.m`
- Signature: `noesyhsqc_ubiquitin_prot()`
- Total lines: 94

## Purpose

1H-1H-15N NOESY-HSQC spectrum of 15N-labelled ubiquitin at 900 MHz with 65 ms mixing time. It is assumed that the protein is not 13C-labelled. Calculation time: hours.

## Physical / mathematical content

- Protein NMR examples. These files specialise liquid-state pulse sequences to labelled biomolecules, exploiting one-bond and two-bond heteronuclear couplings, coherence pathway filtering, selective decoupling, and high-dimensional indirect detection.
- The relaxation model is Redfield-type perturbation theory: fluctuating interactions enter through correlation functions or spectral densities and generate a linear relaxation superoperator.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Implementation structure

- 1H-1H-15N NOESY-HSQC spectrum of 15N-labelled ubiquitin at 900
- MHz with 65 ms mixing time. It is assumed that the protein is
- not 13C-labelled.
- Calculation time: hours.
- Protein data import
- Magnet field
- Tolerances
- Relaxation theory
- Basis set
- Algorithmic options
- Create the spin system structure
- Kill carbons

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `protein()`, `create()`, `kill_spin()`, `strcmp()`, `basis()`, `liquid()`, `apodisation()`, `fftshift()`, `conj()`, `kfigure()`, `plot_3d()`.
