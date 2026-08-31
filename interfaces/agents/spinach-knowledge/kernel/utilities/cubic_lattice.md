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

## Code-derived implementation details

### Comment-guided execution stages

- Lines 29-30: Check consistency; implemented by `grumble(isotope,spacing,n_periods)`.
- Lines 32-33: Generate the isotope list; implemented by `sys.isotopes=cell(1,n_periods^3)`.
- Lines 38-39: Generate coordinates; implemented by `inter.coordinates=cell(n_periods^3,1)`.
- Lines 48-51: Generate lattice translation vectors; implemented by `inter.pbc={spacing*n_periods*[1 0 0], spacing*n_periods*[0 1 0], spacing*n_periods*[0 0 1]}`.

### Control flow inferred from the code

- Line 34: `for` loop over `n=1:n_periods^3`.
- Line 40: `for` loop over `n=1:n_periods`.
- Line 41: `for` loop over `k=1:n_periods`.
- Line 42: `for` loop over `m=1:n_periods`.

### Key state/data transformations

- Lines 33: computes `sys.isotopes` using `sys.isotopes=cell(1,n_periods^3)`.
- Lines 35: computes `sys.isotopes{n}` using `sys.isotopes{n}=isotope`.
- Lines 39: computes `inter.coordinates` using `inter.coordinates=cell(n_periods^3,1)`.
- Lines 49-51: computes `inter.pbc` using `inter.pbc={spacing*n_periods*[1 0 0], spacing*n_periods*[0 1 0], spacing*n_periods*[0 0 1]}`.

### Local helper functions

- Line 56: `grumble()` — `function grumble(isotope,spacing,n_periods)`.
  - Representative operation: `if ~ischar(isotope)`.
  - Representative operation: `error('isotope must be a character string.')`.

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
