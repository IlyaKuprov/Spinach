# examples/relaxation_theory/trosy_fluorine_num.m

- Signature: `trosy_fluorine_num()`

## Purpose

Calculates field-dependent transverse relaxation matrix elements for the `19F` and directly bonded `13C` in the source's 3-fluorotyrosine model. Although the source describes a labelled protein, the constructed spin system contains only this two-spin fragment; it is not a whole-protein simulation.

## Model and quantities

The two spins' coordinates and chemical-shift tensors are selected from `3_fluoro_tyr.log`, using DFT indices 8 for `19F` and 7 for `13C`. The full `sphten-liouv` basis is used without approximation; relaxation is lab-frame Redfield with zero equilibrium and `tau_c = 25e-9` (no unit is annotated for this parameter). The source estimates a calculation time of minutes. At each of 20 fields corresponding to proton Larmor frequencies from 200 to 800 MHz, it obtains the relaxation superoperator and evaluates matrix elements for the single-spin `L+` operators and paired coherence operators. For each nucleus the paired states combine that nucleus's `L+` with plus or minus twice the partner's `Lz` term; the resulting rates represent the broad/narrow TROSY components alongside the single-spin transverse rate.

Two plots show the `19F` and `13C` results. Their horizontal axis is proton Larmor frequency in MHz and their vertical axis is a relaxation matrix element in Hz. This script evaluates relaxation rates directly; it does not simulate or compare an experimental spectrum.

Source: [examples/relaxation_theory/trosy_fluorine_num.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/relaxation_theory/trosy_fluorine_num.m).
