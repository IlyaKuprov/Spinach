# examples/nmr_solids/static_powder_nqi_b.m

- Signature: `static_powder_nqi_b()`

## Purpose

Calculates the static powder 79Br NMR spectrum of potassium bromide. The source notes that at least three quadrupolar tensors are needed to reproduce the experimental line shape, and suggests a distribution of electrostatic environments as a possible explanation. Estimated runtime: seconds.

## Spin system and interactions

The three 79Br sites are at 9.3659 T and each is assigned an isotropic shift of 60.0933 ppm. The quadrupolar matrices are `1e3*diag([13.7569,1.6424,-(13.7569+1.6424)])`, `1e3*diag([4.0779,4.5179,-(4.0779+4.5179)])`, and `1e3*diag([1.5885,0.9449,-(1.5885+0.9449)])`. The calculation uses the `sphten-liouv` basis with IK-0 approximation, projection +1, and inter-level 1; trajectory-level algorithms are disabled.

## Simulation and processing

A 79Br powder acquisition uses the `icos_2ang_163842pts` grid, 1024 points, 100 kHz sweep, 6034.96 Hz receiver offset, and 4096-point zero-fill (frequency axis in Hz, inverted). The initial state is a weighted sum of site-specific 79Br `L+` states with weights 40, 32, and 28; the coil state is the total 79Br `L+` state. The FID receives exponential apodisation with parameter 6 before Fourier transformation and plotting; the displayed vertical range is -10 to 1000.
