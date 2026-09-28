# examples/nmr_nucleic/expt_data/rna_noesy_expt.m

- Signature: `rna_noesy_expt()`

## Purpose

Plots the experimental 2D 1H–1H NOESY spectrum of the Harvard RNA from `rna_noesy_expt.mat` (variable `spec_expt`). This is an experimental-data plotting example, not a simulated pulse-sequence calculation.

## Physical / mathematical content

The spectrum is displayed with its sign reversed. The example sets the field to 17.62 T and defines both frequency axes in ppm, with proton spins in both dimensions.

## Numerical / algorithmic content

The plotting parameters specify offsets of 3473, sweeps of [7500, 7496.252], and zero filling to [1024, 4096] points. The data are passed to `plot_2d` with the source's display and contour settings.

## Implementation structure

The function loads only `spec_expt` from the MAT-file, configures the plot, and calls `plot_2d(spin_system,-spec_expt,...)`. The source credits Shunsuke Imai, Scott Robson, Gerhard Wagner, Zenawi Welderufael, and Ilya Kuprov.
