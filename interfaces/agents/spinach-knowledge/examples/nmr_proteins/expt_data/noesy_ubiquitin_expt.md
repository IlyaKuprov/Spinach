# examples/nmr_proteins/expt_data/noesy_ubiquitin_expt.m

- Signature: `noesy_ubiquitin_expt()`

## Purpose

Plots a precomputed experimental proton spectrum of human ubiquitin loaded from `noesy_ubiquitin_expt.mat`.

## Physical / mathematical content

- Uses proton axes for a two-dimensional spectrum.

## Numerical / algorithmic content

- Loads the spectrum directly and plots it without further spectral processing.

## Implementation structure

- Donghan Lee (Max Planck Institute)
- Ilya Kuprov (University of Southampton)
- Sets the magnetic field to 21.1356 T, offset to 4250 Hz, sweep to 10815 Hz, zero filling to [1024 1024], and axis units to ppm.
- Loads `spectrum` from `noesy_ubiquitin_expt.mat` and plots it with `plot_2d`.
