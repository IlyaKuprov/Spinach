# examples/optimal_control/distortions/distortions_figure_1.m

- Signature: `distortions_figure_1()`

## Purpose

Compare a discretised RF pulse waveform with outputs from simple filter models, showing how pulse-shape distortion affects its in-phase and quadrature components. This is a waveform illustration, not a spin-dynamics calculation or an optimisation.

## Physical / mathematical content

The function converts a pulse table to a two-component waveform and applies single-pole, single-zero, and RLC-style filter models. It does not set a spin system or a spin Hamiltonian. The plotted quadratures are labelled in-phase and quadrature, and the amplitude is expressed as B1 in mT.

## Numerical / algorithmic content

The five-row pulse table spans 0 to 20 microseconds. `fapt2sfo` samples it on a 1000-point grid from -5 to 40 microseconds; the waveform is converted to mT by `1e3*wave/spin('1H')`. For the second-order low-pass comparison, two `spf` applications use `0.9*exp(-1i*0.05)`. For the third-order high-pass comparison, three `szf` applications use `0.1*exp(-1i*0.05)`. The RLC-style comparison applies `spf` twice, with coefficients `0.9*exp(-1i*0.5)` and `0.9*exp(+1i*0.5)`. The function plots each filtered quadrature against the corresponding input component.

## Syntax

Call `distortions_figure_1()` with no arguments. Source: [examples/optimal_control/distortions/distortions_figure_1.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/optimal_control/distortions/distortions_figure_1.m).

## Parameters / inputs

The pulse table, sample grid, conversion, and filter coefficients are defined in the function; it accepts no inputs. Time is displayed in microseconds and B1 in mT.

## Outputs

The function creates three input-versus-output plots for the filter comparisons. It does not return a waveform or a quantitative distortion measure.

## Header notes

The source labels this as Figure 1 from a paper by Rasulov and Kuprov and links [arXiv:2502.02198](https://arxiv.org/abs/2502.02198).
