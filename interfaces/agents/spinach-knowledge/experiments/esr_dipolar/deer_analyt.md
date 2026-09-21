# experiments/esr_dipolar/deer_analyt.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/experiments/esr_dipolar/deer_analyt.m`
- Signature: `deer=deer_analyt(D,J,t)`
- Total lines: 58

## Purpose

Analytical expression for a DEER trace for two spins in the presence of dipolar and exchange coupling. Syntax: deer=deer_analyt(D,J,t)

## Physical / mathematical content

- Dipolar ESR experiment implementations. The pulse logic resolves dipolar couplings by echo modulation, with selective excitation and time-domain accumulation.
- The spin physics includes through-space magnetic dipole-dipole coupling, a rank-2 anisotropic interaction with strong orientation dependence and characteristic secular/non-secular structure.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Parameters / inputs

- D -dipolar coupling, angular frequency units, the
- coefficient in front of (1-3*cos(theta)^2)*Lz*Sz
- in the spin Hamiltonian
- J -exchange coupling, angular frequency units, NMR
- convention (no factor of 2 in front), the coef-
- ficient in front of L*S in the spin Hamiltonian
- t -array of time points, seconds
- Output:
- deer -an array of DEER form factor values of the same
- dimension as t

## Implementation structure

- Analytical expression for a DEER trace for two spins in the
- presence of dipolar and exchange coupling. Syntax:
- deer=deer_analyt(D,J,t)
- D -dipolar coupling, angular frequency units, the
- coefficient in front of (1-3*cos(theta)^2)*Lz*Sz
- in the spin Hamiltonian
- J -exchange coupling, angular frequency units, NMR
- convention (no factor of 2 in front), the coef-
- ficient in front of L*S in the spin Hamiltonian
- t -array of time points, seconds
- Output:
- deer -an array of DEER form factor values of the same

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `fresnelc()`, `fresnels()`, `deer()`, `isscalar()`, `any()`.
