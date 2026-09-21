# kernel/utilities/fpl2phan.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/utilities/fpl2phan.m`
- Signature: `phan=fpl2phan(rho,coil,dims)`
- Total lines: 59

## Purpose

Returns the image painted within the Fokker-Planck vector by the user-specified spin state. Syntax: phan=fpl2phan(rho,coil,dims)

## Physical / mathematical content

- General mathematical and infrastructure utilities. This area contains finite differences, perturbation theory, graph algorithms, spectral densities, tensor algebra, hash/report helpers, and other reusable numerical components.
- The file uses a Fokker-Planck-style enlarged state space in which spatial or orientational coordinates are promoted to extra dimensions and coupled to spin dynamics through differential operators.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Parameters / inputs

- rho -state vector in Fokker-Planck space
- coil -observable state vector in Liouville space
- dims -spatial dimensions of the Fokker-Planck
- problem, a row vector of integers
- Output:
- phan -the image painted by the specified state

## Implementation structure

- Returns the image painted within the Fokker-Planck vector by
- the user-specified spin state. Syntax:
- phan=fpl2phan(rho,coil,dims)
- rho -state vector in Fokker-Planck space
- coil -observable state vector in Liouville space
- dims -spatial dimensions of the Fokker-Planck
- problem, a row vector of integers
- Output:
- phan -the image painted by the specified state
- Check consistency
- Expose the spin dimension
- Compute the observable

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `any()`, `space()`.
