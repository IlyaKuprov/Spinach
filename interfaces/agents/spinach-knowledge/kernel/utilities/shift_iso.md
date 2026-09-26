# kernel/utilities/shift_iso.m

- Signature: `tensors=shift_iso(tensors,spin_numbers,new_iso)`

## Purpose

Replaces the isotropic parts of interaction tensors with user- supplied values. This is useful for correcting DFT calculations, where the anisotropy of the various spin interactions is usual- ly satisfactory, but the isotropic part is not. Syntax: tensors=shift_iso(tensors,spin_numbers,new_iso)

## Physical / mathematical content

- General mathematical and infrastructure utilities. This area contains finite differences, perturbation theory, graph algorithms, spectral densities, tensor algebra, hash/report helpers, and other reusable numerical components.

## Numerical / algorithmic content

## Parameters / inputs

- tensors -a cell array of interaction tensors
- as 3x3 matrices
- spin_numbers -a vector containing the numbers
- of spins in the tensors array that
- should have the isotropic parts
- replaced
- new_iso -a vector containing the new isotro-
- pic parts in the same order as the
- spin numbers listed in spin_numbers

## Outputs

- tensors -a cell array of interaction tensors
- as 3x3 matrices

## Implementation structure

- Replaces the isotropic parts of interaction tensors with user-
- supplied values. This is useful for correcting DFT calculations,
- where the anisotropy of the various spin interactions is usual-
- ly satisfactory, but the isotropic part is not. Syntax:
- tensors=shift_iso(tensors,spin_numbers,new_iso)
- tensors -a cell array of interaction tensors
- as 3x3 matrices
- spin_numbers -a vector containing the numbers
- of spins in the tensors array that
- should have the isotropic parts
- replaced
- new_iso -a vector containing the new isotro-
