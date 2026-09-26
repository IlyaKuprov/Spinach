# kernel/utilities/fpl2rho.m

- Signature: `rho=fpl2rho(rho,dims)`

## Purpose

Returns the average of the spin state vector across the spatial dimensions of the sample. Syntax: rho=fpl2rho(rho,dims)

## Physical / mathematical content

- General mathematical and infrastructure utilities. This area contains finite differences, perturbation theory, graph algorithms, spectral densities, tensor algebra, hash/report helpers, and other reusable numerical components.
- The file uses a Fokker-Planck-style enlarged state space in which spatial or orientational coordinates are promoted to extra dimensions and coupled to spin dynamics through differential operators.

## Numerical / algorithmic content

## Parameters / inputs

- rho -Fokker-Planck state vector
- dims -spatial dimensions of the
- Fokker-Planck problem, a
- vector of positive integers

## Outputs

- rho -Liouville space state vector

## Implementation structure

- Returns the average of the spin state vector across the
- spatial dimensions of the sample. Syntax:
- rho=fpl2rho(rho,dims)
- rho -Fokker-Planck state vector
- dims -spatial dimensions of the
- Fokker-Planck problem, a
- vector of positive integers
- rho -Liouville space state vector
- Check consistency
- Find out the stack size
- Expose the spin dimension (no ND sparse support yet)
- Average over the spatial coordinates
