# kernel/utilities/get_coupling.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/utilities/get_coupling.m`
- Signature: `A=get_coupling(spin_system,n,k)`
- Total lines: 54

## Purpose

Extracts the 3x3 coupling tensor between a pair of spins back from the spin_system data structure. Syntax: A=get_coupling(spin_system,n,k)

## Physical / mathematical content

- General mathematical and infrastructure utilities. This area contains finite differences, perturbation theory, graph algorithms, spectral densities, tensor algebra, hash/report helpers, and other reusable numerical components.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Parameters / inputs

- n,k -indices of the two spins as they appear
- in spin_system.comp.isotopes

## Outputs

- A -3x3 coupling tensor in rad/s

## Implementation structure

- Extracts the 3x3 coupling tensor between a pair of spins back
- from the spin_system data structure. Syntax:
- A=get_coupling(spin_system,n,k)
- n,k -indices of the two spins as they appear
- in spin_system.comp.isotopes
- A -3x3 coupling tensor in rad/s
- Check consistency
- Pull forward and backward coupling
- Fill in empties
- Add up
- Consistency enforcement
- Единственное, что я понимаю в арбузах -это если я по

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `isfield()`, `isscalar()`.
