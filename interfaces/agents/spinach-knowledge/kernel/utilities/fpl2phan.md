# kernel/utilities/fpl2phan.m

- Signature: `phan=fpl2phan(rho,coil,dims)`

## Purpose

Returns the image painted within the Fokker-Planck vector by the user-specified spin state. Syntax: phan=fpl2phan(rho,coil,dims)

## Physical / mathematical content

- Projects the Fokker-Planck vector `rho` onto the coil sensitivity vector with `coil' * rho`, then reshapes the result to dimensions `dims`.

## Numerical / algorithmic content

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
