# examples/nmr_proteins/expt_data/hsqc_ubiquitin_expt.m

- Signature: `hsqc_ubiquitin_expt()`

## Purpose

Processes and plots an experimental HSQC spectrum of human ubiquitin.

## Physical / mathematical content

- Displays a two-dimensional 15N and 1H spectrum.

## Numerical / algorithmic content

- Applies a phase factor and cosine apodisation to positive and negative FIDs, Fourier-transforms and combines them into a States signal, then Fourier-transforms, flips, and plots the spectrum.

## Implementation structure

- Donghan Lee (Max Planck Institute)
- Ilya Kuprov (University of Southampton)
- Loads `hsqc_ubiquitin_expt.mat`, sets the magnetic field to 11.7395 T, zero filling to [1024 1024], sweeps to [2000 4000] Hz, and offsets to [-5870 3753] Hz; plots in ppm.
