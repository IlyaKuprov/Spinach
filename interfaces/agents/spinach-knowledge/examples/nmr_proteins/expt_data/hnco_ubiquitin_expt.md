# examples/nmr_proteins/expt_data/hnco_ubiquitin_expt.m

- Signature: `hnco_ubiquitin_expt()`

## Purpose

Processes and plots an experimental HNCO spectrum of human ubiquitin.

## Physical / mathematical content

- Displays a three-dimensional 15N, 13C, and 1H spectrum.

## Numerical / algorithmic content

- Loads and truncates the FID, applies cosine apodisation, processes F3, F2, and F1 with Fourier transforms, shifts the spectrum, corrects its baseline, and zeros the first 40 points of the third dimension to eliminate the water signal.

## Implementation structure

- Donghan Lee (Max Planck Institute)
- Ilya Kuprov (University of Southampton)
- Loads `hnco_ubiquitin_expt.mat`, sets the magnetic field to 11.7395 T and the spectral axis parameters, then plots the real spectrum in ppm.
