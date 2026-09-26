# kernel/utilities/cubic_lattice.m

- Signature: `[sys,inter]=cubic_lattice(isotope,spacing,n_periods)`

## Purpose

Creates a periodic volume-centered cubic lattice with user- supplied parameters. Syntax: [sys,inter]=cubic_lattice(isotope,spacing,n_periods)

## Physical / mathematical content

- General mathematical and infrastructure utilities. This area contains finite differences, perturbation theory, graph algorithms, spectral densities, tensor algebra, hash/report helpers, and other reusable numerical components.

## Numerical / algorithmic content

## Parameters / inputs

- isotope -character string specifying the isotope,
- for example, '13C'
- spacing -lattice spacing in Angstroms
- n_periods -number of lattice periods in each of the
- three spatial dimensions

## Outputs

- sys, inter -Spinach input data structures with the
- following fields set:
- sys.isotopes, inter.coordinates, inter.pbc

## Implementation structure

- Creates a periodic volume-centered cubic lattice with user-
- supplied parameters. Syntax:
- [sys,inter]=cubic_lattice(isotope,spacing,n_periods)
- isotope -character string specifying the isotope,
- for example, '13C'
- spacing -lattice spacing in Angstroms
- n_periods -number of lattice periods in each of the
- three spatial dimensions
- sys, inter -Spinach input data structures with the
- following fields set:
- sys.isotopes, inter.coordinates, inter.pbc
- Check consistency
