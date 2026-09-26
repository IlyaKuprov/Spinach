# kernel/utilities/phan2fpl.m

- Signature: `rho=phan2fpl(phan,rho)`

## Purpose

Projects a spatial intensity distribution into the Fokker-Planck space, using it as the image painted by the the spin state supp- lied. Syntax: rho=phan2fpl(phan,rho)

## Physical / mathematical content

- General mathematical and infrastructure utilities. This area contains finite differences, perturbation theory, graph algorithms, spectral densities, tensor algebra, hash/report helpers, and other reusable numerical components.
- The file uses a Fokker-Planck-style enlarged state space in which spatial or orientational coordinates are promoted to extra dimensions and coupled to spin dynamics through differential operators.

## Numerical / algorithmic content

## Parameters / inputs

- phan -phantom (the spatial distribution of the
- amplitude of the specified spin state)
- rho -Liouville space state vector

## Outputs

- rho -Fokker-Planck state vector

## Implementation structure

- Projects a spatial intensity distribution into the Fokker-Planck
- space, using it as the image painted by the the spin state supp-
- lied. Syntax:
- rho=phan2fpl(phan,rho)
- phan -phantom (the spatial distribution of the
- amplitude of the specified spin state)
- rho -Liouville space state vector
- rho -Fokker-Planck state vector
- Check consistency
- Stretch the phantom and kron it with the spin state
- Consistency enforcement
- Q: "How many members of a certain demographic group does
