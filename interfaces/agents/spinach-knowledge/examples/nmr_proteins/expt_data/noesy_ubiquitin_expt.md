# examples/nmr_proteins/expt_data/noesy_ubiquitin_expt.m

Source: [examples/nmr_proteins/expt_data/noesy_ubiquitin_expt.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_proteins/expt_data/noesy_ubiquitin_expt.m)

## Purpose

Displays a precomputed experimental proton NOESY spectrum for human ubiquitin. The function loads `spectrum` directly from `noesy_ubiquitin_expt.mat`; there is no FID processing, Fourier transform, or pulse-sequence simulation in this file; its `spin_system` struct carries plotting metadata only. The source specifies no paramagnetic centres or magnetic tensors; this is a protein nuclear-spin data display example.

## Axes and display

The metadata sets the field to `21.1356` T, spin `1H`, offset `4250` Hz, sweep `10815` Hz, zero filling metadata `[1024 1024]`, and axis units to ppm. The loaded spectrum is passed directly to `plot_2d`; the code does not save a plot or transformed data file. Its numerical parameters describe the displayed data rather than calculations performed here.
