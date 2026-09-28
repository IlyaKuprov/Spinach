# examples/esr_liq_pulsed/data_import/orca_import_example.m

- Signature: `orca_import_example()`

## Purpose

Imports ORCA DFT parameters for a methyl-radical ESR simulation. The example attributes its unusual signal-intensity pattern to cross-correlation between the g tensor and hyperfine couplings.

## Physical / mathematical content

The script reads the electron g tensor and three electron–proton hyperfine tensors from the ORCA output, then models the methyl radical with Redfield relaxation and a 0.5 ns correlation time for liquid-state ESR acquisition.

## Numerical / algorithmic content

The acquisition uses 512 points over a 500 MHz sweep; the FID is Fourier transformed with zero filling to 1024 points. The plotting parameters request a derivative spectrum.

## Implementation structure

It parses the ORCA output with `oparse`, transfers the g and hyperfine tensors (converting the couplings to MHz), creates a full sphten-liouv basis, sets electron raising operators for initial and detected states, then runs `liquid` with `@acquire`. It applies the configured apodisation, transforms the FID, and plots the real spectrum.
