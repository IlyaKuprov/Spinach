# kernel/steady.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/steady.m`
- Signature: `rho=steady(spin_system,P,rho,method)`
- Total lines: 220

## Purpose

Steady state under the repeated action by the same dissi- pative evolution propagator. Syntax: rho=steady(spin_system,P,rho,tol,method)

## Physical / mathematical content

- This file belongs to the `kernel` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.
- The optimisation logic is Newton or Newton-like: search directions use first- and second-order local curvature information, usually with regularisation or line-search safeguards.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `size()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Parameters / inputs

- P -propagator, an exponential of the Liouvillian
- that contains a thermalised relaxation super-
- operator (inter.equilibrium='IME' or 'dibari')
- or a product thereof (for example, from a re-
- peating block of a pulse or a pulse sequence)
- rho -optional initial guess for the steady state,
- a good one can significantly accelerate this
- function (leave empty otherwise); the state
- must have unit trace, which in sphten-liouv
- means a first element equal to 1
- method -'newton' (default) for the Newton-Raphson
- steady state solver, 'squaring' for propa-
- gator squaring (much more expensive, but
- unconditionally numerically stable)

## Outputs

- rho -steady state under the repeated applicati-
- on of the propagator P
- Note: available for sphten-liouv and zeeman-liouv formalisms; the
- Newton-Raphson solver pins the first state vector element in
- sphten-liouv and the density matrix trace in zeeman-liouv.

## Implementation structure

- Steady state under the repeated action by the same dissi-
- pative evolution propagator. Syntax:
- rho=steady(spin_system,P,rho,tol,method)
- P -propagator, an exponential of the Liouvillian
- that contains a thermalised relaxation super-
- operator (inter.equilibrium='IME' or 'dibari')
- or a product thereof (for example, from a re-
- peating block of a pulse or a pulse sequence)
- rho -optional initial guess for the steady state,
- a good one can significantly accelerate this
- function (leave empty otherwise); the state
- must have unit trace, which in sphten-liouv

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `exist()`, `rho()`, `speye()`, `complex()`, `grumble()`, `strcmp()`, `clean_up()`, `ismember()`, `ischar()`, `iscolumn()`.
