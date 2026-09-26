# kernel/coherent.m

- Signature: `rho=coherent(spin_system,mode,alpha)`

## Purpose

Coherent state of a bosonic mode. Builds the normalised trunca- tion of the coherent state with the specified amplitude on the specified bosonic mode, with unit operators on all other parti- cles of the system. Syntax: rho=coherent(spin_system,mode,alpha)

## Physical / mathematical content

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

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
