# examples/nmr_nucleic/expt_data/rna_noesy_expt.m

- Signature: `rna_noesy_expt()`
- Source: [examples/nmr_nucleic/expt_data/rna_noesy_expt.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_nucleic/expt_data/rna_noesy_expt.m)

## What this example does

This is a plotting wrapper for an existing experimental 2D 1H-1H NOESY spectrum of the Harvard RNA, as identified in the source comment. It loads the variable `spec_expt` from `rna_noesy_expt.mat` and passes that data to `plot_2d`; it does not calculate the spectrum or define the experiment's pulse program. The MAT-file is an experimental-data input, not a simulation output.

## Display setup

The wrapper sets the field to 17.62 T, identifies both axes with 1H spins, and sets `axis_units='ppm'`. Its remaining axis settings are `offset=3473`, `sweep=[7500.00 7496.252]`, and `zerofill=[1024 4096]`. The source does not annotate units for the offset or sweep values. The plot uses the negative of `spec_expt` and supplies contour/display arguments `20`, `[0.001 0.05 0.001 0.05]`, `2`, `256`, `6`, and `'positive'`; these are plotting settings, not acquisition parameters.

## Scope and limits

No pulse timing, mixing period, gradients, receiver phase/cycling, or experimental processing history is present in this wrapper. It does not inspect or validate the MAT-file contents, report measured peak positions, or compare data with a calculation. The source comment credits Shunsuke Imai, Scott Robson, Gerhard Wagner, Zenawi Welderufael, and Ilya Kuprov.
