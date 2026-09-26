# kernel/utilities/offsetof.m

- Signature: `offs=offsetof(spin_system,idx)`

## Purpose

Returns the isotropic Zeeman offset of the specified spin from the pure magnetogyric ratio frequency in the current magnet. Syntax: offs=offsetof(spin_system,idx)

## Physical / mathematical content

- General mathematical and infrastructure utilities. This area contains finite differences, perturbation theory, graph algorithms, spectral densities, tensor algebra, hash/report helpers, and other reusable numerical components.

## Numerical / algorithmic content

## Parameters / inputs

- idx -index of the spin in sys.isotopes
- array, use idxof() to find index
- by the text label

## Outputs

- offs -offset from the pure magnetogyric
- ratio frequency at the current fi-
- eld, Hz

## Implementation structure

- Returns the isotropic Zeeman offset of the specified spin
- from the pure magnetogyric ratio frequency in the current
- magnet. Syntax:
- offs=offsetof(spin_system,idx)
- idx -index of the spin in sys.isotopes
- array, use idxof() to find index
- by the text label
- offs -offset from the pure magnetogyric
- ratio frequency at the current fi-
- eld, Hz
- Check consistency
- Pull out the Zeeman tensor
