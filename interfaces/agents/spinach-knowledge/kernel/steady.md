# kernel/steady.m

- Signature: `rho=steady(spin_system,P,rho,method)`

## Purpose

Steady state under the repeated action by the same dissi- pative evolution propagator. Syntax: rho=steady(spin_system,P,rho,tol,method)

## Physical / mathematical content

- The optimisation logic is Newton or Newton-like: search directions use first- and second-order local curvature information, usually with regularisation or line-search safeguards.

## Numerical / algorithmic content

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
