# examples/nmr_proteins/noesyhsqc_ubiquitin_deut.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_proteins/noesyhsqc_ubiquitin_deut.m`
- Signature: `noesyhsqc_ubiquitin_deut()`
- Total lines: 101

## Purpose

1H-1H-15N NOESY-HSQC spectrum of 15N-labelled ubiquitin at 900 MHz with 90 ms mixing time. It is assumed that the protein is not 13C-labelled. Specific positions are deuterated, and deu- terium nuclei are simulated explicitly as spin-1 particles. Calculation time: a week on 32 cores, needs 512 GB of RAM.

## Physical / mathematical content

- Protein NMR examples. These files specialise liquid-state pulse sequences to labelled biomolecules, exploiting one-bond and two-bond heteronuclear couplings, coherence pathway filtering, selective decoupling, and high-dimensional indirect detection.
- The relaxation model is Redfield-type perturbation theory: fluctuating interactions enter through correlation functions or spectral densities and generate a linear relaxation superoperator.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Implementation structure

- 1H-1H-15N NOESY-HSQC spectrum of 15N-labelled ubiquitin at 900
- MHz with 90 ms mixing time. It is assumed that the protein is
- not 13C-labelled. Specific positions are deuterated, and deu-
- terium nuclei are simulated explicitly as spin-1 particles.
- Calculation time: a week on 32 cores, needs 512 GB of RAM.
- Protein data import
- Magnet field
- Tolerances
- Relaxation theory
- Basis set
- Algorithmic options
- Create the spin system structure

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `protein()`, `num2cell()`, `strcmp()`, `create()`, `kill_spin()`, `basis()`, `liquid()`, `apodisation()`, `fftshift()`, `conj()`, `kfigure()`, `plot_3d()`.
