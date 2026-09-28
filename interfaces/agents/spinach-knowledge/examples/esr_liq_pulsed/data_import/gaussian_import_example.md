# examples/esr_liq_pulsed/data_import/gaussian_import_example.m

- Signature: `gaussian_import_example()`

## Purpose

Imports Gaussian DFT parameters for a methyl-radical ESR simulation. The example attributes its unusual signal-intensity pattern to cross-correlation between the g tensor and hyperfine couplings.

## Physical / mathematical content

The imported system contains an electron and a proton, with the isotope mapping supplied to `g2spinach`. The example uses a Redfield relaxation model with a 0.5 ns correlation time and simulates liquid-state ESR acquisition.

## Numerical / algorithmic content

The acquisition uses 512 points over a 500 MHz sweep; the FID is Fourier transformed with zero filling to 1024 points. The plotting parameters request a derivative spectrum.

## Implementation structure

It parses the Gaussian output through `gparse` and `g2spinach`, creates a full sphten-liouv basis, sets electron raising operators for initial and detected states, then runs `liquid` with `@acquire`. It applies the configured apodisation, transforms the FID, and plots the real spectrum.
