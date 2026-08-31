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

## Code-derived implementation details

### Comment-guided execution stages

- Lines 46-47: Check consistency; implemented by `grumble(group_name)`.
- Lines 49-50: Do the necessary; implemented by `switch group_name`.
- Lines 320-321: Complain and bomb out; implemented by `error(['permutation group ' group_name ' is not available.'])`.
- Lines 325-326: Get the list of group elements; implemented by `group.elements=vertcat(group.class{:})`.
- Lines 328-329: Transform class-wise character list into element-wise list; implemented by `group.characters=zeros(length(group.class_sizes),sum(group.class_sizes))`.

### Control flow inferred from the code

- Line 50: dispatches on `group_name`; cases `'S2'`, `'S3'`, `'S4'`, `'S4A'`, `'S5'`, `'S6'`, `'S6A'`.
- Line 330: `for` loop over `n=1:length(group.class_sizes)`.
- Line 331: `for` loop over `k=1:group.class_sizes(n)`.

### Key state/data transformations

- Lines 54: computes `group.name` using `group.name='S2 permutation group (Abelian)'`.
- Lines 55: computes `group.order` using `group.order=2`.
- Lines 56: computes `group.nclasses` using `group.nclasses=2`.
- Lines 57: computes `group.class_sizes` using `group.class_sizes=[1 1]`.
- Lines 58: computes `group.class{1}` using `group.class{1}=[1 2]`.
- Lines 59: computes `group.class{2}` using `group.class{2}=[2 1]`.
- Lines 60: computes `group.n_irreps` using `group.n_irreps=2`.
- Lines 61: computes `group.irrep_dims` using `group.irrep_dims=[1 1]`.
- Lines 62: computes `group.class_characters` using `group.class_characters=[1 1`.
- Lines 73: computes `group.class{3}` using `group.class{3}=[1 3 2; 3 2 1; 2 1 3]`.
- Lines 89: computes `group.class{4}` using `group.class{4}=[1 2 4 3; 1 3 2 4; 1 4 3 2; 2 1 3 4; 4 2 3 1; 3 2 1 4]`.
- Lines 90: computes `group.class{5}` using `group.class{5}=[2 1 4 3; 3 4 1 2; 4 3 2 1]`.
- Lines 138: computes `group.class{6}` using `group.class{6}=[2 1 4 5 3; 2 1 5 3 4; 2 3 1 5 4; 2 4 5 1 3; 2 5 4 3 1`.
- Lines 142: computes `group.class{7}` using `group.class{7}=[2 3 4 5 1; 2 3 5 1 4; 2 4 1 5 3; 2 4 5 3 1; 2 5 1 3 4`.
- Lines 211: computes `group.class{8}` using `group.class{8}=[2 1 4 3 6 5; 2 1 5 6 3 4; 2 1 6 5 4 3; 3 4 1 2 6 5; 3 5 1 6 2 4; 3 6 1 5 4 2; 4 3 2 1 6 5; 4 5 6 1 2 3; 4 6 5 1 3 2; 5 3 2 6 1 4`.
- Lines 213: computes `group.class{9}` using `group.class{9}=[2 3 1 5 6 4; 2 3 1 6 4 5; 2 4 5 1 6 3; 2 4 6 1 3 5; 2 5 4 6 1 3; 2 5 6 3 1 4; 2 6 4 5 3 1; 2 6 5 3 4 1; 3 1 2 5 6 4; 3 1 2 6 4 5`.
- Lines 217: computes `group.class{10}` using `group.class{10}=[2 1 4 5 6 3; 2 1 4 6 3 5; 2 1 5 3 6 4; 2 1 5 6 4 3; 2 1 6 3 4 5; 2 1 6 5 3 4; 2 3 4 1 6 5; 2 3 5 6 1 4; 2 3 6 5 4 1; 2 4 1 3 6 5`.
- Lines 226: computes `group.class{11}` using `group.class{11}=[2 3 4 5 6 1; 2 3 4 6 1 5; 2 3 5 1 6 4; 2 3 5 6 4 1; 2 3 6 1 4 5; 2 3 6 5 1 4; 2 4 1 5 6 3; 2 4 1 6 3 5; 2 4 5 3 6 1; 2 4 5 6 1 3`.

### Local helper functions

- Line 339: `grumble()` — `function grumble(group_name)`. It's so wonderful to see a great, new, crucial idea which is not mine!
  - Representative operation: `if ~ischar(group_name)`.
  - Representative operation: `error('group_name must be a character string.')`.

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
