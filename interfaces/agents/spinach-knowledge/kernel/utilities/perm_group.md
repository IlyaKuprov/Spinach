# kernel/utilities/perm_group.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/utilities/perm_group.m`
- Signature: `group=perm_group(group_name)`
- Total lines: 349

## Purpose

Permutation group database. Returns complete data for permutation groups. The following group names are available: S2, S3, S4, S4A, S5, S6, S6A, S8A. The options ending in A are the largest real-valued-character Abelian subgroups. Syntax: group=perm_group(group_name)

## Physical / mathematical content

- General mathematical and infrastructure utilities. This area contains finite differences, perturbation theory, graph algorithms, spectral densities, tensor algebra, hash/report helpers, and other reusable numerical components.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Parameters / inputs

- group_name -character string specifying the name
- group, e.g. 'S5'

## Outputs

- group.name -long name of the group
- group.order -number of elements in the group
- group.nclasses -number of classes in the group
- group.class_sizes -a row vector giving number of elements
- in each class
- group.class -a cell array of matrices giving the ele-
- ments belonging to each class. The ele-
- ments are given as row vectors of permu-
- tation strings stacked vertically into
- a matrix.
- group.n_irreps -number of irreducible representations in
- the group
- group.irrep_dims -a row vector giving dimensions of irre-
- ducible representation
- group.class_characters -a matrix of characters for each irredu-
- cible representation (in rows) of each
- class (in columns).

## Implementation structure

- Permutation group database. Returns complete data for permutation
- groups. The following group names are available: S2, S3, S4,
- S4A, S5, S6, S6A, S8A. The options ending in A are the
- largest real-valued-character Abelian subgroups. Syntax:
- group=perm_group(group_name)
- group_name -character string specifying the name
- group, e.g. 'S5'
- group.name -long name of the group
- group.order -number of elements in the group
- group.nclasses -number of classes in the group
- group.class_sizes -a row vector giving number of elements
- in each class

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `group()`, `vertcat()`, `ischar()`.
