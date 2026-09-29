# examples/dnp_liq/ccdnp/rates_si_sys_b.m

- MATLAB implementation: [examples/dnp_liq/ccdnp/rates_si_sys_b.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/dnp_liq/ccdnp/rates_si_sys_b.m)

## What it computes

Builds a Redfield relaxation superoperator for the system-B, one-proton/two-electron cross-correlated DNP model and prints selected proton/electron self- and cross-relaxation projections. It is a rate-inspection example, not a field sweep or a time-domain DNP calculation.

## Running and inputs

Call the no-argument MATLAB function `rates_si_sys_b()` with Spinach on the MATLAB path. All system values are hard-coded; the script reads no external inputs. It calls Spinach to create the system, make a full `sphten-liouv` basis (no approximation), construct the Redfield superoperator, and form the operator states used in the printout. The source comments estimate seconds of calculation time.

## System B parameterisation

- Isotopes: `{'1H','E','E'}`; field assignment `sys.magnet=14.1` (unit not stated in the source). Proton Zeeman values are `[0 10 20]` with Euler triple `[0 0 0]`.
- Electron 1 Zeeman values are `[1.977800 1.977600 1.977600]` with Euler triple `[-0.180 0.017 0.194]`. Electron 2 values are `[2.006800 2.003800 2.003800]` with Euler triple `[0.632 0.783 1.086]`. The scalar electron-electron coupling is assigned `5e6`.
- Coordinates are proton `[0 0 0]`, electron 1 `[6.000 0.030 0.317]`, and electron 2 `[-6.000 -0.038 0.535]`. The source does not identify coordinate units.
- Relaxation uses `inter.relaxation={'redfield'}`, `inter.equilibrium='zero'`, and `inter.rlx_keep='labframe'`. It sets `inter.temperature=298`, `inter.tau_c={100e-12}`, and `sys.tols.rlx_integration=1e-10`; no units are stated for these literals.

## Output and interpretation

After `R=relaxation(spin_system)`, the script prints normalised longitudinal self-relaxation projections for proton and both electrons, electron-to-proton projections, four transverse mixed-state projections (`E1p ± 2*E1pE2z` and `E2p ± 2*E1zE2p`), and proton-containing longitudinal cross terms, including the terms labelled for NzE1z, NzE2z, and NzE1zE2z to Nz. It also prints mixed transverse-coherence projections labelled E1p to NzE1p, E1p to NzE1pE2z, and E1pE2z to NzE1pE2z. A distinctive point in this system-B setup is the combination of electron-1 values near 1.9778 with electron-2 values near 2.0068/2.0038 and separate Euler triples; it is not the near-matched electron tensor pair used by system A. These printed contractions are selected elements/projections, not a complete relaxation matrix. No files or figures are produced.

## Reference and caveat

The source cites https://doi.org/10.1016/j.jmr.2021.106940. The source does not annotate units for the magnet, Zeeman, coupling, coordinate, temperature, correlation-time, or tolerance values, so the settings are quoted as written rather than assigned inferred units.
