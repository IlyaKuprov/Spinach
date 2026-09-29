# kernel/overloads/@ttclass/revert.m

- Signature: `tt=revert(tt)`

## Action

For each train column and each core, permute the core axes as `[4,2,3,1]`, then reverse the core sequence. For a core shaped `[r_left,m,n,r_right]`, its permuted shape is `[r_right,m,n,r_left]`: the two bond axes exchange sides, while the row and column mode axes remain in place. The core count is unchanged; each core's row and column sizes remain in their original axes, while the ordered mode sequence and bond-rank sequence reverse with the cores. The coefficient field is not assigned or scaled by this method.

This is the source-described bit-revert permutation. It transforms the TT cores and leaves the result in TT form; it does not assemble the represented matrix. No conjugation or row/column transpose is performed; the method adds no explicit input guard.

## Input and output

- `tt` — tensor-train operator.
- `tt` — the same train collection with reversed core order and bond direction.

## Source

- [MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/overloads/@ttclass/revert.m)
- [Spin Dynamics Wiki](https://spindynamics.org/wiki/index.php?title=ttclass/revert.m)

D. Savostyanov.
