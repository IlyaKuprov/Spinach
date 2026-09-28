# kernel/utilities/fpl2rho.m

- Signature: `rho=fpl2rho(rho,dims)`

## Purpose

Returns the average of the spin state vector across the spatial dimensions of the sample. Syntax: rho=fpl2rho(rho,dims)

## Physical / mathematical content

- Reshapes each column of the Fokker-Planck state vector to separate spin and spatial dimensions, averages across the spatial dimension, and squeezes the result.

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
