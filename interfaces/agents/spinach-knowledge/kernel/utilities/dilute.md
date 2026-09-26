# kernel/utilities/dilute.m

- Signature: `subsystems=dilute(spin_system,isotope,tuples)`

## Purpose

Splits the spin system into several independent subsystems, each containing only one instance of a user specified isotope or iso- tope tuple that is deemed "dilute". All spin system data is upda- ted accordingly. Basis set information, if found, is destroyed.

## Physical / mathematical content

- General mathematical and infrastructure utilities. This area contains finite differences, perturbation theory, graph algorithms, spectral densities, tensor algebra, hash/report helpers, and other reusable numerical components.

## Numerical / algorithmic content

## Syntax

```matlab
subsystems=dilute(spin_system,isotope,tuples)
```

## Parameters / inputs

- spin_system -spin system object obtained
- from create.m function
- isotope -isotope specification string,
- for example '13C'
- tuples -1 returns all subsystems with
- a single instance of the isoto-
- pe, 2 returns all subsystems
- where two instances of the iso-
- tope are present, etc.

## Outputs

- subsystems -a cell array of spin system
- objects, with each cell cor-
- responding to one of the new-
- ly formed isotopomers.

## Implementation structure

- Splits the spin system into several independent subsystems, each
- containing only one instance of a user specified isotope or iso-
- tope tuple that is deemed "dilute". All spin system data is upda-
- ted accordingly. Basis set information, if found, is destroyed.
- subsystems=dilute(spin_system,isotope,tuples)
- spin_system -spin system object obtained
- from create.m function
- isotope -isotope specification string,
- for example '13C'
- tuples -1 returns all subsystems with
- a single instance of the isoto-
- pe, 2 returns all subsystems
