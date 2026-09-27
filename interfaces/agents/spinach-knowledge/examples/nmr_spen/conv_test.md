# examples/nmr_spen/conv_test.m

- Signature: `conv_test()`

## Purpose

Tests convergence of spatial spin dynamics for a pulsed-field-gradient spin echo by varying spatial grid size and the finite-difference derivative stencil. It fits the simulated attenuation to the Stejskal–Tanner relation to estimate diffusion and plots the estimate and runtime. The source reports minutes on an NVIDIA Tesla A100 and substantially longer on CPU.

## Model and scan

The model is one 1H spin at 11.7426 T, with no relaxation, in a 15 mm sample. The reference diffusion coefficient is 18×10⁻¹⁰ m²/s; the small gradient duration is 2 ms and the diffusion delay is 50 ms. A Gaussian initial spatial profile and zero velocity field are used. Fifty grid sizes span 1000–10000 points, three periodic finite-difference stencils (3, 5, and 7 points) are tested, and 20 gradient amplitudes span 0–0.5 T/m.

For each grid/stencil pair the PFG spin echo is simulated, intensities are normalised to the zero-gradient signal, and the fitted diffusion estimate is compared with the reference value. The figures show the diffusion estimate and elapsed time versus grid size, and attenuation curves across gradient amplitudes.
