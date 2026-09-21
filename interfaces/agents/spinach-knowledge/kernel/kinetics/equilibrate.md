# kernel/kinetics/equilibrate.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/kinetics/equilibrate.m`
- Signature: `c=equilibrate(K,c0)`
- Total lines: 95

## Purpose

Equilibrates linear chemical kinetics and returns a vector of equilibrium concentrations. Syntax: c=equilibrate(K,c0)

## Physical / mathematical content

- This file belongs to the `kernel` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `size()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

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

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `scomponents()`, `logical()`, `vertcat()`, `cond()`, `eps()`, `all()`, `iscolumn()`, `any()`.
