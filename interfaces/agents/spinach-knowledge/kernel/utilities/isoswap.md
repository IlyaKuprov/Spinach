# kernel/utilities/isoswap.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/utilities/isoswap.m`
- Signature: `[sys,inter]=isoswap(sys,inter,spins,new_iso)`
- Total lines: 116

## Purpose

Makes isotope replacements in the input structures. All interactions are automatically scaled as necessary. Syntax: [sys,inter]=isoswap(sys,inter,spins,new_iso)

## Physical / mathematical content

- General mathematical and infrastructure utilities. This area contains finite differences, perturbation theory, graph algorithms, spectral densities, tensor algebra, hash/report helpers, and other reusable numerical components.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Parameters / inputs

- sys, inter -Spinach input data structures
- spins -a vector of integers specifying
- spin numbers
- new_iso -character string specifying
- the new isotope
- Output:
- sys, inter -Spinach input data structures
- Note: quadratic and higher order couplings are wiped and a warning
- is printed -those are not transferable.

## Implementation structure

- Makes isotope replacements in the input structures. All interactions
- are automatically scaled as necessary. Syntax:
- [sys,inter]=isoswap(sys,inter,spins,new_iso)
- sys, inter -Spinach input data structures
- spins -a vector of integers specifying
- spin numbers
- new_iso -character string specifying
- the new isotope
- Output:
- Note: quadratic and higher order couplings are wiped and a warning
- is printed -those are not transferable.
- Grumbler missing

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `isfield()`, `spins()`, `int2str()`, `wiped()`, `spin()`, `setdiff()`, `isstruct()`, `isvector()`, `any()`, `ischar()`.
