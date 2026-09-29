# examples/nmr_solids/static_powder_dip.m

- Signature: `static_powder_dip()`
- Source: [examples/nmr_solids/static_powder_dip.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_solids/static_powder_dip.m)

## Purpose

A static-powder pulse-acquire calculation for a dipolar-coupled two-spin system; the source records a calculation time of seconds.

## Spin model and powder average

The pair comprises two 1H spins with scalar Zeeman entries 5.0 and -2.0 and coordinates `[0 0 0]` and `[0 3.9 0.1]` (coordinate units are not stated). The field parameter is 14.1, with no unit given. The static powder average uses `rep_2ang_6400pts_sph`; the example specifies no rotor or gradient sequence.

## Acquisition and processing

The 1H channel starts from and detects `L+`, with no decoupling. Acquisition uses offset 1000, sweep 12000, 128 points, and zero filling to 512. The display axis is configured in ppm; units for offset and sweep are not stated. The FID is exponentially apodised with parameter 6, Fourier transformed with the specified zero fill, and plotted as its real part.
