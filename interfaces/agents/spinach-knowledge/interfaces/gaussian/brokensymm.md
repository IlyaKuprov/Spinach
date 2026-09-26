# interfaces/gaussian/brokensymm.m

- Signature: `J=brokensymm(props_sing,props_trip)`

## Purpose

Exchange coupling estimation from a pair of DFT logs using Yamaguchi equation. The notation is: H=-2J*(Sa.Sb)

## Physical / mathematical content

- Gaussian interfaces. These parse quantum-chemistry output into spin Hamiltonian ingredients such as hyperfine, shielding, or exchange parameters.
- The relevant state manifold is the singlet/triplet decomposition, where permutation symmetry controls selection rules, relaxation susceptibility, and convertibility to ordinary magnetisation.

## Numerical / algorithmic content

## Syntax

```matlab
J=brokensymm(props_sing,props_trip)
```

## Parameters / inputs

- props_sing -the output of gparse for the singlet
- state of the biradical
- props_trip -the output of gparse for the triplet
- state of the biradical

## Outputs

- J -an order-of-magnitude (really rough)
- estimate of exchange coupling, Hz

## Implementation structure

- Exchange coupling estimation from a pair of DFT logs using
- Yamaguchi equation. The notation is:
- H=-2J*(Sa.Sb)
- J=brokensymm(props_sing,props_trip)
- props_sing -the output of gparse for the singlet
- state of the biradical
- props_trip -the output of gparse for the triplet
- J -an order-of-magnitude (really rough)
- estimate of exchange coupling, Hz
- Check consistency
- Eq 6 in https://doi.org/10.1063/1.5144696
- Convert from Hartree to Hz
