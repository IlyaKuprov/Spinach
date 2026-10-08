# examples/nmr_liquids/hetcor_strychnine.m

- MATLAB implementation: [examples/nmr_liquids/hetcor_strychnine.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_liquids/hetcor_strychnine.m)

- Signature: `hetcor_strychnine()`

## Purpose and spin system

Magnitude-mode liquid-state HETCOR for strychnine at natural 13C abundance; the source comment estimates minutes of calculation time. The wrapper requests `strychnine({'1H','13C'})`, so its active network is proton and carbon spin sites (the helper defines 22 `1H` and 21 `13C` sites; its two `15N` sites are not selected). It uses `dilute(...,'13C')` to form the natural-abundance carbon isotopomers. The helper supplies shifts and scalar couplings from Berger and Braun, except the one-bond C18-H18b coupling (http://dx.doi.org/10.1016/j.jmr.2014.02.003), and coordinates for a major conformer (http://dx.doi.org/10.1039/C0CC04114A). The helper has no 13C-13C couplings, which it identifies as a limitation to natural-abundance 13C simulations.

## Acquisition and simulation

The wrapper sets `sys.magnet=5.9`, enables `zte` and `greedy` with proximity cutoff 4.0, and uses `sphten-liouv` / `IK-2`, scalar-coupling connectivity, and proximity level 1. Its sequence settings are `J=140` Hz, sweeps `[3000 10000]` Hz, offsets `[1000 4000]` (units not stated here), `[256 256]` points, and `[512 512]` zero-fill points; `spins={'1H','13C'}` makes F1 proton and F2 carbon, with `1H` decoupling and ppm axes. The separate `experiments/nmr_liquids/hetcor.m` pulse program starts from proton longitudinal magnetisation, detects on carbon, and implements the magnitude-mode HETCOR transfer with fixed delays `1/(2J)` and `1/(3J)`; the wrapper supplies parameters rather than spelling out those pulse operations. The program receives the Hamiltonian, relaxation, and kinetics superoperators from the `liquid(...,'nmr')` context; this wrapper sets no relaxation rates or model.

The wrapper simulates each diluted isotopomer in `parfor`, applies cosine apodisation in both dimensions, sums shifted 2D Fourier transforms, and plots the absolute spectrum in positive mode. The HETCOR pulse-program reference is https://doi.org/10.1016/0022-2364(81)90272-9.
