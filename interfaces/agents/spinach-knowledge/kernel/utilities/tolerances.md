# kernel/utilities/tolerances.m

- Signature: `[spin_system,sys]=tolerances(spin_system,sys)`

## Purpose

Tolerances and fundamental constants. Sets various accuracy cut-offs, constants and tolerances used by Spinach kernel. Syntax: spin_system=tolerances(spin_system,sys)

## Physical / mathematical content

- General mathematical and infrastructure utilities. This area contains finite differences, perturbation theory, graph algorithms, spectral densities, tensor algebra, hash/report helpers, and other reusable numerical components.
- The spin physics includes through-space magnetic dipole-dipole coupling, a rank-2 anisotropic interaction with strong orientation dependence and characteristic secular/non-secular structure.

## Numerical / algorithmic content

- A Krylov-subspace or Arnoldi construction is used to avoid forming or exponentiating very large dense propagators directly.

## Parameters / inputs

- spin_system -Spinach system description object
- sys -system specification object described
- in the input preparation section of
- the manual

## Outputs

- spin_system -updated system description object
- sys -system specification structure with
- the tolerance substructure parsed out
- Notes: direct calls and modifications to this function are discouraged:
- the accuracy settings should be modified by setting the sys.tols
- structure, see the input preparation manual.

## Implementation structure

- Tolerances and fundamental constants. Sets various accuracy cut-offs,
- constants and tolerances used by Spinach kernel. Syntax:
- spin_system=tolerances(spin_system,sys)
- spin_system - Spinach system description object
- sys - system specification object described
- in the input preparation section of
- the manual
- spin_system - updated system description object
- sys - system specification structure with
- the tolerance substructure parsed out
- the accuracy settings should be modified by setting the sys.tols
- structure, see the input preparation manual.
