# kernel/utilities/overwound.m

- Signature: `overwound(rho,spc_dim,spn_dim)`

## Purpose

Checks whether a Fokker–Planck state vector contains spatial frequencies that its grid is dangerously close to misrepresenting because it has too few points.

## Parameters

- `rho` — Fokker–Planck state vector or a bookshelf stack of state vectors. It must be a numeric array with `prod([spn_dim spc_dim])` rows; its columns form the stack.
- `spc_dim` — spatial dimensions `[X Y Z]`, specified as three positive integers.
- `spn_dim` — spin dimension, specified as a positive integer.

## Diagnostics

The function checks input consistency, then examines each spatial dimension whose point count exceeds one. An examined dimension with fewer than 10 points produces an error requesting a higher point count in that dimension. Otherwise, the function applies `fft` and `fftshift` along that dimension, sums the absolute Fourier amplitudes over the remaining dimensions and stack, and plots the result. Each plot has a spatial-frequency axis relative to the Nyquist limit, running from `-1` to `1`, and a population-density axis in arbitrary units. Output consists of figures and diagnostic messages to the console; the function has no return value.

For spatial dynamics such as diffusion and flow with finite-difference derivative operators, set the spatial grid point count to several times the minimum Nyquist value.

## Reference

- [Spin Dynamics: overwound.m](https://spindynamics.org/wiki/index.php?title=overwound.m)
- Contact: ilya.kuprov@weizmann.ac.il