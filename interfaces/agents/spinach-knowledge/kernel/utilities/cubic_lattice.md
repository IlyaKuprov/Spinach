# kernel/utilities/cubic_lattice.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/utilities/cubic_lattice.m`
- Signature: `[sys,inter]=cubic_lattice(isotope,spacing,n_periods)`
- Total lines: 75

## Purpose

Creates a periodic volume-centered cubic lattice with user- supplied parameters. Syntax: [sys,inter]=cubic_lattice(isotope,spacing,n_periods)

## Physical / mathematical content

- General mathematical and infrastructure utilities. This area contains finite differences, perturbation theory, graph algorithms, spectral densities, tensor algebra, hash/report helpers, and other reusable numerical components.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

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

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `sub2ind()`, `ischar()`.
