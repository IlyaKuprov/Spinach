# kernel/utilities/magpump.m

- Signature: `R=magpump(spin_system,R,rho,rate)`

## Purpose

Adds phenomenological pumping terms to the relaxation superoperator to enable approximate simulation of CIDNP, PHIP and DNP type effects.

## Physical / mathematical content

- Adds phenomenological pumping terms to a relaxation superoperator for approximate CIDNP, PHIP, and DNP simulations.

## Numerical / algorithmic content

- Updates the first column of `R` by adding `rate*rho`; the supported formalism is `sphten-liouv`, and the pumped state must have zero unit-state component.

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

- Checks the matrix, state vector, finite scalar rate, formalism, and unit-state component, then performs the first-column update. Repeated calls can add multiple pumped states.
