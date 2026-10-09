# etc/molecules/dac_reaction.m

- MATLAB implementation: [etc/molecules/dac_reaction.m](https://github.com/IlyaKuprov/Spinach/blob/main/etc/molecules/dac_reaction.m)

**Call:** `[sys, inter, bas, kin] = dac_reaction()`

**Inputs:** none. Run with Spinach available and the function's companion Gaussian log files reachable beside the source: the routine locates its own directory and reads `cyclopentadiene.log`, `acrylonitrile.log`, `norbornene_endo.log`, and `norbornene_exo.log`. Set `sys.magnet` in the caller before a field-dependent simulation; this builder does not set it.

## What it builds

The example merges four proton spin systems and appends a spin-free solvent pool in this order: cyclopentadiene (spin indices 1–6), acrylonitrile (7–9), endo-norbornene carbonitrile (10–18), exo-norbornene carbonitrile (19–27), and acetonitrile (empty spin list, solvent). The header calls the first reactant “pentadiene”, but the implementation explicitly loads `cyclopentadiene.log`; the code-level identity is therefore cyclopentadiene. For the four Gaussian-derived species, `gparse` and `g2spinach` import calculated geometry and chemical-shift anisotropy, after which the code replaces isotropic shifts and scalar couplings with tabulated values. Some coupling signs are explicitly noted as missing in the source.

Acetonitrile is a concentration-only pool; solvent proton shifts, couplings, and relaxation are not represented. The reacting-flow initial state did not excite solvent, and this representation excludes its NMR signal explicitly. Each of the five chemical parts has a unit reference concentration; callers supply their actual initial concentrations or spatial unit-coordinate fields.

## Outputs and settings

- `sys`, `inter`: merged Spinach system and interaction structures.
- `bas`: `formalism='sphten-liouv'`, `approximation='IK-2'`, `connectivity='scalar_couplings'`, `prox_level=1` on spin-bearing parts; the solvent uses `approximation='none'`.
- `kin`: two reactant-to-product matching records, both using reactant parts [1 2]. The first targets part 3 (endo) with matches [1 12; 2 17; 3 18; 4 16; 5 10; 6 11; 7 14; 8 15; 9 13]. The second targets part 4 (exo) with [1 21; 2 26; 3 27; 4 25; 5 19; 6 20; 7 23; 8 24; 9 22]. Each record has additive closure and an initial rate of zero; callers must set the two physical rates in `inter.chem.reactions` before `create`. The fourth output is an independent value copy of that pair: modifying `kin` does not update `inter`. Either set rates directly in `inter.chem.reactions`, or assign `inter.chem.reactions=kin` after modifying the returned records.

The spin-bearing molecules use secular Redfield relaxation and zero equilibrium. Correlation times, in part order, are 5, 20, 50, 50, and 1 ps. The former T1/T2 term acted only on solvent protons; it is absent from the spin-free solvent model. Part 3 and the first reaction are endo; part 4 and the second are exo.

**Source reference:** [Spinach Wiki: dac_reaction.m](https://spindynamics.org/wiki/index.php?title=dac_reaction.m).
