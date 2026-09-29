# examples/optimal_control/distortions/distortions_figure_2.m

- Signature: `distortions_figure_2()`

## Purpose

Illustrate amplitude compression by comparing the `E1000B` pulse with outputs from two amplifier models. The function plots nonlinear waveform transformations; it does not perform spin dynamics or pulse optimisation.

## Physical / mathematical content

The input is a two-channel representation containing the generated waveform and a zero-valued second component. The tanh and root compression models act on this representation with a nominal saturation level of 3e4 rad/s. The plotted ordinate is nutation frequency in rad/s.

## Numerical / algorithmic content

`vg_pulse('E1000B',1000,0.001)` supplies the waveform. A 1000-point time axis covers 0 to 0.001 seconds (0 to 1 ms). The function applies `amp_tanh` with parameter `3e4` and `amp_root` with parameters `3e4` and `10`, then plots the input and each model's first output component. Dashed reference levels are drawn at +3e4 and -3e4 rad/s.

## Syntax

Call `distortions_figure_2()` with no arguments. Source: [examples/optimal_control/distortions/distortions_figure_2.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/optimal_control/distortions/distortions_figure_2.m).

## Parameters / inputs

The pulse name, generator arguments, sampling, and amplifier parameters are fixed in the function; there are no function inputs. The horizontal axis is in milliseconds and the vertical axis is in rad/s.

## Outputs

The function creates a plot of the input, tanh-model output, root-model output (legend label `s=10`), and the two saturation levels. It does not return a waveform or a measured amplifier response.

## Header notes

The source labels this as Figure 2 from a paper by Rasulov and Kuprov and links [arXiv:2502.02198](https://arxiv.org/abs/2502.02198).
