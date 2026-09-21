# kernel/utilities/rlx_scalar.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/utilities/rlx_scalar.m`
- Signature: `R=rlx_scalar(spin_system,H0,H1,tau_c_array)`
- Total lines: 90

## Purpose

Scalar relaxation superoperator using Redfield theory. Syntax: R=rlx_scalar(spin_system,H0,H1,tau_c_array)

## Physical / mathematical content

- General mathematical and infrastructure utilities. This area contains finite differences, perturbation theory, graph algorithms, spectral densities, tensor algebra, hash/report helpers, and other reusable numerical components.
- The relaxation model is Redfield-type perturbation theory: fluctuating interactions enter through correlation functions or spectral densities and generate a linear relaxation superoperator.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Parameters / inputs

- H0 -background Hamiltonian
- H1 -the stochastically modulated interaction operator
- multiplied by its root mean square modulation depth
- tau_c_array -a cell array of the following format:
- {[weight_a,tau_a],[weight_b,tau_b],...}
- giving weights of the exponential components
- of the correlation function and the associa-
- ted correlation times, e.g. {[1.0,1e-12]}

## Outputs

- R -relaxation superoperator, a negative definite matrix
- Note: if H1(t) has a non-zero time or ensemble average value,
- that average must be subtracted out and placed into H0

## Implementation structure

- Scalar relaxation superoperator using Redfield theory. Syntax:
- R=rlx_scalar(spin_system,H0,H1,tau_c_array)
- H0 -background Hamiltonian
- H1 -the stochastically modulated interaction operator
- multiplied by its root mean square modulation depth
- tau_c_array -a cell array of the following format:
- {[weight_a,tau_a],[weight_b,tau_b],...}
- giving weights of the exponential components
- of the correlation function and the associa-
- ted correlation times, e.g. {[1.0,1e-12]}
- R -relaxation superoperator, a negative definite matrix
- Note: if H1(t) has a non-zero time or ensemble average value,

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `clean_up()`, `expmint()`, `speye()`, `ismatrix()`, `ishermitian()`, `all()`, `iscell()`.
