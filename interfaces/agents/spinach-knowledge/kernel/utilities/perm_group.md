# kernel/utilities/perm_group.m

**Source:** <https://github.com/IlyaKuprov/Spinach/blob/main/kernel/utilities/perm_group.m>

## Purpose

`perm_group` is a permutation group database that returns complete data for a requested permutation group. The available group names are `S2`, `S3`, `S4`, `S4A`, `S5`, `S6`, `S6A`, and `S8A`. The options ending in `A` are the largest real-valued-character Abelian subgroups.

## Behaviour

The function is called as `group=perm_group(group_name)`. It first validates the input through an internal consistency check (`grumble`), which errors with `'group_name must be a character string.'` if the argument is not a character string. A `switch` statement then selects the group data; an unrecognised name triggers the error `'permutation group ' group_name ' is not available.'`.

After the switch, the function assembles `group.elements` by vertically concatenating all class matrices (`vertcat(group.class{:})`), and expands the class-wise character table into an element-wise character matrix `group.characters`, where each column repeats the character value of the corresponding class for every element in that class.

Stored data per group includes:

- `S2`: order 2, 2 classes of sizes `[1 1]`, 2 irreps of dimensions `[1 1]`.
- `S3`: order 6, 3 classes of sizes `[1 2 3]`, 3 irreps of dimensions `[1 1 2]`.
- `S4`: order 24, 5 classes of sizes `[1 6 8 6 3]`, 5 irreps of dimensions `[1 1 2 3 3]`.
- `S4A` (largest Abelian subgroup of S4): order 4, 4 classes of sizes `[1 1 1 1]`, 4 irreps of dimensions `[1 1 1 1]`.
- `S5`: order 120, 7 classes of sizes `[1 10 20 15 30 20 24]`, 7 irreps of dimensions `[1 1 4 4 6 5 5]`.
- `S6`: order 720, 11 classes of sizes `[1 15 40 45 90 120 144 15 40 90 120]`, 11 irreps of dimensions `[1 1 5 5 10 10 9 9 5 5 16]`.
- `S6A` (largest real-valued-character Abelian subgroup of S6): order 8, 8 classes of sizes `[1 1 1 1 1 1 1 1]`, 8 irreps of dimensions `[1 1 1 1 1 1 1 1]`.
- `S8A` (largest real-valued-character Abelian subgroup of S8): order 16, 16 classes of sizes `[1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1]`, 16 irreps of dimensions `[1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1]`.

## Inputs and outputs

**Input:**

- `group_name` — character string specifying the name of the group, e.g. `'S5'`.

**Outputs:**

- `group.name` — long name of the group.
- `group.order` — number of elements in the group.
- `group.nclasses` — number of classes in the group.
- `group.class_sizes` — row vector giving the number of elements in each class.
- `group.class` — cell array of matrices giving the elements belonging to each class; elements are given as row vectors of permutation strings stacked vertically into a matrix.
- `group.n_irreps` — number of irreducible representations in the group.
- `group.irrep_dims` — row vector giving dimensions of irreducible representations.
- `group.class_characters` — matrix of characters for each irreducible representation (in rows) of each class (in columns).
- `group.elements` — all group elements, obtained by vertical concatenation of the class matrices.
- `group.characters` — element-wise character matrix, with each class's character column repeated according to the class sizes.

## References

- Spinach Wiki page for `perm_group.m`: <https://spindynamics.org/wiki/index.php?title=perm_group.m>
- Spinach Dynamics homepage: <https://spindynamics.org/>
