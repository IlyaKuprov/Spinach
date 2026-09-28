# kernel/utilities/perm_group.m

- Signature: `group=perm_group(group_name)`

## Purpose

Returns the stored permutation-group data for the supported groups S2, S3, S4, S4A, S5, S6, S6A, and S8A. The A-suffixed options are the largest real-valued-character Abelian subgroups.

## Parameters

- `group_name` — character string naming the group, for example `'S5'`.

## Outputs

- `group.name` — long group name.
- `group.order` — number of elements.
- `group.nclasses` — number of conjugacy classes.
- `group.class_sizes` — row vector giving the number of elements in each class.
- `group.class` — cell array of matrices containing the permutation rows for each class.
- `group.n_irreps` — number of irreducible representations.
- `group.irrep_dims` — row vector of irreducible-representation dimensions.
- `group.class_characters` — character matrix, with irreducible representations in rows and classes in columns.
- `group.elements` — class matrices concatenated into one matrix.
- `group.characters` — class characters expanded to the individual elements in each class.

An unsupported group name raises an error.

## Implementation structure

The function selects the stored data for the requested group, then forms the full element list from the class matrices and expands the class-character columns across the elements in their respective classes.

## Reference

- <https://spindynamics.org/wiki/index.php?title=perm_group.m>
