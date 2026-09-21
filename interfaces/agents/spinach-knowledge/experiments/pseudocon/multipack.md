# experiments/pseudocon/multipack.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/experiments/pseudocon/multipack.m`
- Signature: `Ilm=multipack(ranks,moments)`
- Total lines: 66

## Purpose

Packs multipole moments from a linear stream into a cell array that is arranged by ranks. Syntax: Ilm=multipack(ranks,moments)

## Physical / mathematical content

- Paramagnetic-pseudocontact inference routines. The mathematics includes inverse problems, tensor parameterisation, interpolation, and regularisation.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Parameters / inputs

- ranks -a vector of spherical ranks present, e.g. [0 1 2]
- moments -a vector of multipole moments for each rank,
- arranged in a linear stream

## Outputs

- Ilm -a cell array of vector corresponding to the multipole
- moments defined in http://dx.doi.org/10.1039/c6cp05437d
- for L=0, one element
- for L=1, three elements
- for L=2, five elements
- et cetera.

## Implementation structure

- Packs multipole moments from a linear stream into a cell array
- that is arranged by ranks. Syntax:
- Ilm=multipack(ranks,moments)
- ranks -a vector of spherical ranks present, e.g. [0 1 2]
- moments -a vector of multipole moments for each rank,
- arranged in a linear stream
- Ilm -a cell array of vector corresponding to the multipole
- moments defined in http://dx.doi.org/10.1039/c6cp05437d
- for L=0, one element
- for L=1, three elements
- for L=2, five elements
- et cetera.

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `moments()`, `ranks()`, `isrow()`, `any()`.
