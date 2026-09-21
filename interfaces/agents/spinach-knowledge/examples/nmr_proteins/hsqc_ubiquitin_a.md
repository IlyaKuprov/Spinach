# examples/nmr_proteins/hsqc_ubiquitin_a.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_proteins/hsqc_ubiquitin_a.m`
- Signature: `hsqc_ubiquitin_a()`
- Total lines: 72

## Purpose

1H-15N HSQC of human ubiquitin, decoupling applied in both dimensions. Calculation time: hours, faster with a Tesla A100 GPU. Zenawi Welderufael Luke Edwards Ilya Kuprov

## Physical / mathematical content

- Protein NMR examples. These files specialise liquid-state pulse sequences to labelled biomolecules, exploiting one-bond and two-bond heteronuclear couplings, coherence pathway filtering, selective decoupling, and high-dimensional indirect detection.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Implementation structure

- 1H-15N HSQC of human ubiquitin, decoupling applied
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
- Spinach housekeeping

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `protein()`, `create()`, `basis()`, `liquid()`, `apodisation()`, `fftshift()`, `conj()`, `kfigure()`, `scale_figure()`, `plot_2d()`.
