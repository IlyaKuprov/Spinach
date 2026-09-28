# examples/optimal_control/distortions/restrans_test.m

- Signature: `restrans_test()`

## Purpose

Tests the resonator transform by sending a square pulse through a simple resonator model and plotting the time-domain response.

## Method

Generates a random 51-sample, two-component pulse and zeros the first and last three samples. It then simulates the response for ¹H and ¹⁵N at 14.1 T, using a 1 μs slice duration, Q factor 50, and piecewise-constant modelling with 100 points. The circuit resonance frequencies are 2π × 600 MHz for ¹H and 2π × 60 MHz for ¹⁵N.
