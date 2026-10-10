# examples/dnp_liq/ccdnp/rates_si_sys_a.m

- MATLAB implementation: [examples/dnp_liq/ccdnp/rates_si_sys_a.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/dnp_liq/ccdnp/rates_si_sys_a.m)

## What it computes

Constructs the Redfield relaxation superoperator for the system-A case of cross-correlated liquid-state DNP, then prints selected relaxation projections for proton and electron operators. The model contains one proton coupled to two exchange-connected electrons; unlike a DNP enhancement simulation, this function reports rate contractions and does not propagate a driven steady state.

## Running and inputs

Run the no-argument MATLAB function `rates_si_sys_a()` with Spinach on the MATLAB path. The script defines the spin system and all settings inline, creates its basis and relaxation superoperator, and needs no external data file or user argument. Its source comment estimates the calculation at seconds. It uses `sphten-liouv` with no basis approximation and Redfield relaxation.

## System A parameterisation

- Isotopes: `{'1H','E','E'}`; magnet assignment `sys.magnet=14.1` T. Proton Zeeman shifts `[0 10 20]` are in ppm, with Euler angles `[0 0 0]` in radians.
- The dimensionless electron g-tensor values are `[1.977873 1.977798 1.977792]` with zero Euler angles for electron 1, and `[1.977919 1.978000 1.978000]` with Euler triple `[-0.590 -0.100 0.490]` radians for electron 2. The scalar exchange assignment is `6.2e6` Hz.
- The file explicitly labels its coordinate set “Coordinates for anisotropic HF”: proton `[0 0 0]`, electron 1 `[7.0300 0.0187 0.9820]`, electron 2 `[-7.0300 0.2051 -1.0001]`. Spinach interprets these coordinates in ångström.
- The relaxation configuration is Redfield with equilibrium mode `zero` and `rlx_keep='labframe'`. It sets temperature `298` K, correlation time `100e-12` s, and dimensionless integration tolerance `1e-10`.

## What it prints

The script computes `R=relaxation(spin_system)` and prints projections for normalised longitudinal proton/electron operators, the electron transverse combinations `E1p ± 2*E1pE2z` and `E2p ± 2*E1zE2p`, and proton/electron longitudinal cross terms, including the terms labelled for NzE1z, NzE2z, and NzE1zE2z to Nz. It also prints mixed transverse-coherence projections labelled E1p to NzE1p, E1p to NzE1pE2z, and E1pE2z to NzE1pE2z. A useful source-specific distinction is that the coordinate comment identifies these coordinates as inputs for anisotropic hyperfine interactions; they are not merely display geometry. The function emits text to the MATLAB command window and does not save rates or make plots.

## Reference and caveat

The source cites https://doi.org/10.1016/j.jmr.2021.106940. The physical units above follow the Spinach parameter contracts, including where this driver omits explicit unit annotations.
