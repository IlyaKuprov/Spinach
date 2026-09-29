# examples/nmr_solids/static_powder_nqi_b.m

- Signature: `static_powder_nqi_b()`
- Source: [examples/nmr_solids/static_powder_nqi_b.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_solids/static_powder_nqi_b.m)

## Purpose and model

Calculates the static powder 79Br NMR spectrum of potassium bromide; the source estimates seconds. Its comment says at least three quadrupolar tensors are needed to reproduce the experimental line shape and suggests, tentatively, a distribution of electrostatic environments as a possible reason.

The three 79Br spins share the isotropic shift 60.0933 ppm at field parameter 9.3659. Their separate quadrupolar matrices are `1e3*diag([13.7569,1.6424,-(13.7569+1.6424)])`, `1e3*diag([4.0779,4.5179,-(4.0779+4.5179)])`, and `1e3*diag([1.5885,0.9449,-(1.5885+0.9449)])`. The basis is `sphten-liouv` with IK-0 approximation, projection +1, and inter-level 1; trajectory-level algorithms are disabled. The static powder average uses `icos_2ang_163842pts`; no rotor or gradient parameters are set.

## Acquisition and processing

The acquisition uses 79Br, sweep 1e5 Hz, receiver offset 6034.96 Hz, 1024 points, and 4096-point zero-fill; the axis is in Hz and inverted. The initial state combines site-specific `L+` states with weights 40, 32, and 28; the coil state is the total 79Br `L+`. After exponential apodisation with parameter 6, the real Fourier spectrum is plotted with vertical limits -10 to 1000.
