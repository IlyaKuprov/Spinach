# examples/relaxation_theory/trosy_fluorine_sym.m

- Signature: `trosy_fluorine_sym()`

## Purpose

An analytical dipole-dipole/chemical-shift-anisotropy calculation of field-dependent TROSY relaxation rates for the `19F` and bonded `13C` in the source's 3-fluorotyrosine model. Coordinates and shift tensors are taken from `3_fluoro_tyr.log`; the source passes them to `rlx_dd_csa` with correlation-time parameter `25e-9` (no unit is annotated for this parameter). The source estimates a calculation time of seconds.

## Rates and plotted mechanisms

The calculation samples 20 fields corresponding to proton Larmor frequencies from 200 to 800 MHz. It plots broad and narrow TROSY line rates together with the total transverse rate for fluorine and carbon; the axes are labelled proton Larmor frequency in MHz and relaxation matrix element in Hz. A third plot decomposes the carbon TROSY rate into dipole-dipole, CSA, and DD-CSA cross-correlation contributions. The source plots the cross term as `-abs(c_tro_xc)`, making the interference contribution visible alongside the two positive mechanism terms.

These are calculated relaxation rates and mechanism contributions from the analytical function, not measured line widths or an experimental spectrum. This example calls `rlx_dd_csa`; it does not set up a Bloch-Redfield or stochastic-Liouville spin-system calculation in this file.

Source: [examples/relaxation_theory/trosy_fluorine_sym.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/relaxation_theory/trosy_fluorine_sym.m).
