# kernel/coherent.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/coherent.m`
- Signature: `rho=coherent(spin_system,mode,alpha)`
- Total lines: 95

## Purpose

Coherent state of a bosonic mode. Builds the normalised trunca- tion of the coherent state with the specified amplitude on the specified bosonic mode, with unit operators on all other parti- cles of the system. Syntax: rho=coherent(spin_system,mode,alpha)

## Physical / mathematical content

- This file belongs to the `kernel` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Parameters / inputs

- mode -index of a bosonic mode in sys.isotopes
- alpha -coherent state amplitude, a complex scalar

## Outputs

- rho -coherent state density matrix (zeeman-hilb)
- or its vectorisation (zeeman-liouv)
- Note: the Fock space truncation of the mode chops the tail of
- the Poisson distribution; the state is renormalised after
- the truncation and the lost weight is reported.

## Implementation structure

- Coherent state of a bosonic mode. Builds the normalised trunca-
- tion of the coherent state with the specified amplitude on the
- specified bosonic mode, with unit operators on all other parti-
- cles of the system. Syntax:
- rho=coherent(spin_system,mode,alpha)
- mode -index of a bosonic mode in sys.isotopes
- alpha -coherent state amplitude, a complex scalar
- rho -coherent state density matrix (zeeman-hilb)
- or its vectorisation (zeeman-liouv)
- Note: the Fock space truncation of the mode chops the tail of
- the Poisson distribution; the state is renormalised after
- the truncation and the lost weight is reported.

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `factorial()`, `report()`, `num2str()`, `conj()`, `speye()`, `rho()`, `isfield()`, `basis()`, `isscalar()`, `ismember()`.
