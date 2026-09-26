# kernel/kinetics/equilibrate.m

- Signature: `c=equilibrate(K,c0)`

## Purpose

Equilibrates linear chemical kinetics and returns a vector of equilibrium concentrations. Syntax: c=equilibrate(K,c0)

## Physical / mathematical content

## Numerical / algorithmic content

## Parameters / inputs

- K -reaction rate matrix corresponding to
- dc/dt=K*c, where c is the concentration
- vector
- c0 -vector of initial concentrations

## Outputs

- c -vector of equilibrium concentrations

## Implementation structure

- Equilibrates linear chemical kinetics and returns a vector of
- equilibrium concentrations. Syntax:
- c=equilibrate(K,c0)
- K -reaction rate matrix corresponding to
- dc/dt=K*c, where c is the concentration
- vector
- c0 -vector of initial concentrations
- c -vector of equilibrium concentrations
- Check consistency
- Shortcut for zero concentrations
- Recursive calls for independent reactions
- Assemble the steady state system
