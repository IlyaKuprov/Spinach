# examples/optimal_control/distortions/distortions_figure_1.m

- Signature: `distortions_figure_1()`

## Purpose

Figure 1 from the paper by Rasulov and Kuprov:

## Physical / mathematical content

- A discretised pulse sequence is converted to a waveform in mT. Cascaded single-pole and single-zero filters, plus an RLC filter, are applied to compare the input and filtered in-phase and quadrature components.

## Numerical / algorithmic content

## Implementation structure

- Figure 1 from the paper by Rasulov and Kuprov:
- Pulse sequence and its discretisation
- Convert waveform to mT
- Apply a cascade of two single-pole filters
- Original vs second-order low-pass filter
- Apply a cascade of three single-zero filters
- Original vs third-order high-pass filter
- Apply an RLC filter
- Original vs RLC filter
