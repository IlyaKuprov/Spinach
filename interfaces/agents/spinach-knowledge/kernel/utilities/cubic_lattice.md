# kernel/utilities/cubic_lattice.m

- Signature: `[sys,inter]=cubic_lattice(isotope,spacing,n_periods)`

## Purpose

Creates a periodic volume-centered cubic lattice with user-supplied parameters.

## Physical / mathematical content

The lattice contains `n_periods^3` sites, with coordinates on a cubic grid separated by `spacing` Angstroms. Periodic translation vectors extend `spacing*n_periods` along each of the three coordinate axes.

## Numerical / algorithmic content

The function validates its inputs, assigns `isotope` to every site, generates the site coordinates, and sets the three periodic translation vectors.

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

- `grumble` checks that `isotope` is a character string, `spacing` is a positive finite real scalar, and `n_periods` is a positive finite integer scalar.
- The function fills `sys.isotopes` with `n_periods^3` copies of `isotope`.
- Nested loops populate `inter.coordinates` with `spacing*[(n-1) (k-1) (m-1)]` for `n`, `k`, and `m` from `1` to `n_periods`.
- `inter.pbc` contains the three axis-aligned translation vectors.