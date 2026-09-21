# etc/molecules/dac_reaction.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/etc/molecules/dac_reaction.m`
- Signature: `[sys,inter,bas,kin]=dac_reaction()`
- Total lines: 194

## Purpose

Example Diels-Alder cycloaddition reaction settings: pentadiene (reactant), acrylonitrile (reactant), exo-norbornene (product), endo-norbornene (product), and acetonitrile (solvent). Atom co- ordinates and chemical shift anisotropies are pulled from a DFT calculation; isotropic chemical shifts and J-couplings are ex- perimental (some J-coupling signs are missing). Syntax: [sys,inter,bas,kin]=dac_reaction()

## Physical / mathematical content

- This file belongs to the `etc` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.
- The relaxation model is Redfield-type perturbation theory: fluctuating interactions enter through correlation functions or spectral densities and generate a linear relaxation superoperator.

## Numerical / algorithmic content

## Parameters / inputs

- none

## Outputs

- sys, inter, bas -Spinach input data structures, remember
- to specify the field in sys.magnet
- kin -matching tables for which nuclei go where in which
- of the two chemical reactions

## Implementation structure

- Example Diels-Alder cycloaddition reaction settings: pentadiene
- (reactant), acrylonitrile (reactant), exo-norbornene (product),
- endo-norbornene (product), and acetonitrile (solvent). Atom co-
- ordinates and chemical shift anisotropies are pulled from a DFT
- calculation; isotropic chemical shifts and J-couplings are ex-
- perimental (some J-coupling signs are missing). Syntax:
- [sys,inter,bas,kin]=dac_reaction()
- none
- sys, inter, bas -Spinach input data structures, remember
- to specify the field in sys.magnet
- kin -matching tables for which nuclei go where in which
- of the two chemical reactions

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `mfilename()`, `own_folder()`, `gparse()`, `g2spinach()`, `shift_iso()`, `idxof()`, `num2cell()`, `merge_inp()`.
