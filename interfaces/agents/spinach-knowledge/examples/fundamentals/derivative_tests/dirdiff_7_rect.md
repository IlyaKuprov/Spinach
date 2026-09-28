# examples/fundamentals/derivative_tests/dirdiff_7_rect.m

- Signature: `dirdiff_7_rect()`

## Purpose

Checks directional derivatives returned by phase-modulated GRAPE for a rectangular pulse and an ensemble of chained waveform distortions.

## Setup

The test runs with sphten-liouv, zeeman-liouv, and zeeman-hilb. For each formalism it builds the system with dirdiff_test_system and sets a 13C control channel, channel map [1;1], drift H, controls Lx and Ly, initial states Sx, Sy, Sz, and targets −Sz, Sy, Sx. The power levels are 2*pi*linspace(50e3,70e3,10); GRAPE uses lbfgs, up to 1000 iterations, and the rectangle integrator.

There are five pulse intervals, each 12.8e-6 s, with unit amplitudes. Two distortion chains are applied: firf(w,[0.9 0.1i]) → spf(w,0.2) → szf(w,0.2) → amp_root(w,2*pi*20e3,4), and the reverse ordering szf(w,0.2) → spf(w,0.2) → amp_root(w,2*pi*20e3,4) → firf(w,[0.9 0.1i]).

## Derivative check

For a random five-sample phase vector randn(1,5)/3, the code obtains the analytical gradient from grape_phase and estimates derivatives with centered finite differences using h=sqrt(eps('double')). It checks samples 1, 5, and 3 (left edge, right edge, and midpoint). Each must satisfy abs(grad_anl-grad_num)/abs(grad_num)<1e-6; a failed comparison raises an error identifying the formalism and sample position.
