# etc/molecules/dac_reaction.m

- Signature: `[sys,inter,bas,kin]=dac_reaction()`

## Purpose

Example Diels-Alder cycloaddition reaction settings: pentadiene (reactant), acrylonitrile (reactant), exo-norbornene (product), endo-norbornene (product), and acetonitrile (solvent). Atom co- ordinates and chemical shift anisotropies are pulled from a DFT calculation; isotropic chemical shifts and J-couplings are ex- perimental (some J-coupling signs are missing). Syntax: [sys,inter,bas,kin]=dac_reaction()

## Physical / mathematical content

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
