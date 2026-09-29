# examples/nmr_liquids/hmbc_cyprinol.m

- MATLAB implementation: [examples/nmr_liquids/hmbc_cyprinol.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_liquids/hmbc_cyprinol.m)

- Signature: `hmbc_cyprinol()`

## Purpose and spin system

Liquid-state HMBC for cyprinol at natural 13C abundance; the source comment estimates seconds of calculation time. The wrapper obtains `sys`, `inter`, and `bas` from `cyprinol()`; that helper defines 42 proton sites and 27 carbon sites (labelled `1H` and `13C`) and says shifts and J-couplings come from http://dx.doi.org/10.1002/mrc.4782, with missing values estimated by the source's own stated method of "tossing a twenty-sided coin." The wrapper applies `dilute(...,'13C')` for the natural-abundance isotopomers.

## Acquisition and simulation

The wrapper sets `sys.magnet=11.7`, enables `greedy` with proximity cutoff 4.0, and uses `sphten-liouv` / `IK-1`, scalar-coupling connectivity, proximity level 1, and inter-level 3. It sets `J=150` Hz, `delta_b=60e-3` s, sweeps `[12000 2500]` Hz, offsets `[5000 1250]` (units not stated here), `[128 128]` points, and `[512 512]` zero-fill points. `spins={'13C','1H'}` assigns carbon to F1 and proton to F2; axes are in ppm.

The wrapper calls `liquid(...,@hmbc,parameters,'nmr')`; the separate `experiments/nmr_liquids/hmbc.m` pulse program excites and detects protons, uses carbon pulses and J-set delays, and evolves carbon in F1 with proton decoupling. Hamiltonian, relaxation, and kinetics superoperators come from the `liquid` context; the wrapper specifies no relaxation rates or model. It simulates diluted isotopomers in `parfor`, applies cosine apodisation in both dimensions, sums shifted 2D Fourier transforms, and plots the absolute spectrum in positive mode. HMBC pulse-program references: https://doi.org/10.1021/ja00268a061 and https://doi.org/10.1016/0022-2364(88)90172-2.
