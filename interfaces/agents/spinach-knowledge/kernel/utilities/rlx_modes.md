# kernel/utilities/rlx_modes.m

- Signature: `R=rlx_modes(spin_system)`

## Purpose

Bosonic mode dissipation superoperator. Builds thermalised GKSL dissipators for the amplitude damping and the pure dephasing of the bosonic modes declared in inter.modes, using the amplitude damping rates and the pure dephasing rates ingested by create.m and the Bose-Einstein thermal occupation numbers computed from the physical mode frequencies, meaning the sum of the declared carrier and the declared frequency wher

## Physical / mathematical content

- General mathematical and infrastructure utilities. This area contains finite differences, perturbation theory, graph algorithms, spectral densities, tensor algebra, hash/report helpers, and other reusable numerical components.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Parameters / inputs

- spin_system -Spinach spin system description object
- with bosonic mode information present

## Outputs

- R -bosonic mode dissipation superoperator
- Note: the dissipators are kappa*(1+nbar)*D[a], kappa*nbar*D[c],
- and 2*gamma_phi*D[n], where D[x] is the GKSL dissipator of
- the operator x, built from ladder operators truncated to
- the level count of each mode. The Spinach convention of
- zero temperature meaning the high-temperature limit is not
- applicable to bosonic modes: zero temperature here produces
- zero thermal occupation numbers.

## Implementation structure

- Bosonic mode dissipation superoperator. Builds thermalised GKSL
- dissipators for the amplitude damping and the pure dephasing of
- the bosonic modes declared in inter.modes, using the amplitude
- damping rates and the pure dephasing rates ingested by create.m
- and the Bose-Einstein thermal occupation numbers computed from
- the physical mode frequencies, meaning the sum of the declared
- carrier and the declared frequency where inter.modes.carriers
- is present, and the system temperature. Syntax:
- R=rlx_modes(spin_system)
- spin_system -Spinach spin system description object
- with bosonic mode information present
- R -bosonic mode dissipation superoperator
