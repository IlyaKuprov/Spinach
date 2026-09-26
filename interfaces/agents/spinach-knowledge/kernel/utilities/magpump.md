# kernel/utilities/magpump.m

- Signature: `R=magpump(spin_system,R,rho,rate)`

## Purpose

Adds phenomenological pumping terms to the relaxation superoperator to enable approximate simulation of CIDNP, PHIP and DNP type effects.

## Physical / mathematical content

- General mathematical and infrastructure utilities. This area contains finite differences, perturbation theory, graph algorithms, spectral densities, tensor algebra, hash/report helpers, and other reusable numerical components.

## Numerical / algorithmic content

## Syntax

```matlab
R=magpump(spin_system,R,rho,rate)
```

## Parameters / inputs

- R -relaxation superoperator, from relaxation()
- rho -the state to be pumped, from state()
- rate -pumping rate, Hz

## Outputs

- R -modified relaxation superoperator
- Note: for the pumping to work correctly, the unit state population
- (first element) in the state vector that R will be acting on
- must be set to 1.
- Note: this function is only available in sphten-liouv formalism, and
- may be called repeatedly if multiple states are pumped.

## Implementation structure

- Adds phenomenological pumping terms to the relaxation superoperator
- to enable approximate simulation of CIDNP, PHIP and DNP type effects.
- R=magpump(spin_system,R,rho,rate)
- R -relaxation superoperator, from relaxation()
- rho -the state to be pumped, from state()
- rate -pumping rate, Hz
- R -modified relaxation superoperator
- Note: for the pumping to work correctly, the unit state population
- (first element) in the state vector that R will be acting on
- must be set to 1.
- Note: this function is only available in sphten-liouv formalism, and
- may be called repeatedly if multiple states are pumped.
