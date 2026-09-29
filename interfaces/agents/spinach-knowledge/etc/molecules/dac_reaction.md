# etc/molecules/dac_reaction.m

- MATLAB implementation: [etc/molecules/dac_reaction.m](https://github.com/IlyaKuprov/Spinach/blob/main/etc/molecules/dac_reaction.m)

**Call:** `[sys, inter, bas, kin] = dac_reaction()`

**Inputs:** none. Run with Spinach available and the function's companion Gaussian log files reachable beside the source: the routine locates its own directory and reads `cyclopentadiene.log`, `acrylonitrile.log`, `norbornene_endo.log`, and `norbornene_exo.log`. Set `sys.magnet` in the caller before a field-dependent simulation; this builder does not set it.

## What it builds

The example assembles five proton spin systems and merges them in this order: cyclopentadiene (spin indices 1–6), acrylonitrile (7–9), endo-norbornene carbonitrile (10–18), exo-norbornene carbonitrile (19–27), and acetonitrile (28–30, solvent). The header calls the first reactant “pentadiene”, but the implementation explicitly loads `cyclopentadiene.log`; the code-level identity is therefore cyclopentadiene. For the four Gaussian-derived species, `gparse` and `g2spinach` import calculated geometry and chemical-shift anisotropy, after which the code replaces isotropic shifts and scalar couplings with tabulated values. Some coupling signs are explicitly noted as missing in the source.

Acetonitrile is instead specified directly: three proton shifts are set to 2.0, its coordinates are empty, and its scalar-coupling matrix is zero. Each compound is a separate chemical part with a concentration entry of 1; these are example relative entries, not calibrated experimental concentrations.

## Outputs and settings

- `sys`, `inter`: merged Spinach system and interaction structures.
- `bas`: `formalism='sphten-liouv'`, `approximation='IK-2'`, `connectivity='scalar_couplings'`, `prox_level=1`.
- `kin`: two reactant-to-product matching records, both using reactant parts [1 2]. The first targets part 3 (endo) with matches [1 12; 2 17; 3 18; 4 16; 5 10; 6 11; 7 14; 8 15; 9 13]. The second targets part 4 (exo) with [1 21; 2 26; 3 27; 4 25; 5 19; 6 20; 7 23; 8 24; 9 22]. These are atom/spin correspondences, not reaction rate constants.

The relaxation configuration is Redfield with T1/T2 rates retained in the secular approximation and zero equilibrium; correlation times, in part order, are 5, 20, 50, 50, and 1 ps. R1/R2 entries are initialised to zero, with entries 28–30 set to 0.5 for solvent; the source does not state rate units here. The comment beside the second product assignment mistakenly calls part 4 “endo”; its actual definition above is exo, so the index and matching table identify the exo product.

**Source reference:** [Spinach Wiki: dac_reaction.m](https://spindynamics.org/wiki/index.php?title=dac_reaction.m).
