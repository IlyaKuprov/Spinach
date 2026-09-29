# examples/nmr_solids/static_powder_csa.m

- Signature: `static_powder_csa()`
- Source: [examples/nmr_solids/static_powder_csa.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_solids/static_powder_csa.m)

## Purpose

A two-spin static CSA powder pattern; the source records a calculation time of seconds.

## Spin model and powder average

The system contains two 1H spins with distinct anisotropic Zeeman eigenvalue inputs, `[-2 -2 4]-5` and `[-1 -3 4]+5`; both Euler-angle triples are `[0 0 0]`. The source supplies no units for these entries. The field parameter is 14.1, also without an explicit unit. The static powder calculation uses `rep_2ang_6400pts_sph`; no rotor or gradient sequence is configured.

## Acquisition and processing

The 1H channel starts from and detects `L+`, with no decoupling. Acquisition uses offset 0, sweep 15000, 256 points, and zero filling to 512. The display axis is configured in ppm; units for offset and sweep are not stated. The FID is exponentially apodised with parameter 6, Fourier transformed with the specified zero fill, and plotted as its real part.
