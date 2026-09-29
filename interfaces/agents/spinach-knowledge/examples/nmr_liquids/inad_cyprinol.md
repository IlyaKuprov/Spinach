# examples/nmr_liquids/inad_cyprinol.m

- MATLAB implementation: [examples/nmr_liquids/inad_cyprinol.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_liquids/inad_cyprinol.m)

- Signature: `inad_cyprinol()`

## Purpose

Build a 1D INADEQUATE spectrum from naturally abundant 13C pairs in cyprinol. The source comment estimates minutes for the calculation.

## Spin system and method

`cyprinol()` supplies 42 `1H` and 27 `13C` sites with scalar-coupling and isotropic-shift data; its source says unreported values are estimated. The script forms two-13C isotopomers using `dilute(spin_system,'13C',2)`, computes the pair coupling as `trace(get_coupling(...))/3`, and simulates only pairs satisfying `abs(J)>2*pi*1.0` (the source comment describes this as stronger than 1 Hz). It uses a sparse Liouville-space basis (`sphten-liouv`, IK-1), scalar-coupling connectivity, proximity level 1 and interaction level 4. Greedy setup uses proximity cutoff 4.0; the field is 11.7 T.

## INADEQUATE acquisition

Only `13C` is active, with `1H` decoupling. The wrapper specifies a working `J=50 Hz`, sweep 10000 Hz, offset 5000 (the sources do not state its unit), 4096 points and 8192 zero-fill points; the axis is in ppm and inverted. Its source comment describes selection of double-quantum coherence from coupled 13C pairs and conversion back for detection.

## Processing and scope

Each accepted pair is simulated with `liquid(...,@inadequate,...,'nmr')`; the FID receives exponential apodisation with the source value 6, is Fourier-transformed and added to the 1D spectrum. The real spectrum is plotted. No explicit relaxation parameters are set in the wrapper. It chooses the eligible isotope pairs, method/basis and acquisition/processing settings; pulse-train details are in `experiments/nmr_liquids/inadequate.m`.

## Sources

Cyprinol shift/coupling source: http://dx.doi.org/10.1002/mrc.4782
INADEQUATE sequence source: https://doi.org/10.1021/ja00534a056
