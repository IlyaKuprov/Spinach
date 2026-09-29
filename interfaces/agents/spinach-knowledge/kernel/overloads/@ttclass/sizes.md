# kernel/overloads/@ttclass/sizes.m

- Signature: `modesizes=sizes(tt)`

## Action

For each core in the first train, the method records the second and third core dimensions. It returns an `ncores`-by-2 array whose row `k` is `[size(tt.cores{k,1},2), size(tt.cores{k,1},3)]`, i.e. that core's physical row and column dimensions. It does not return bond ranks or aggregate dimensions across cores.

This is a metadata query: it leaves the TT cores and ranks unchanged and does not materialise the represented matrix. It applies no conjugation or transpose and adds no explicit input guard.

## Input and output

- `tt` — tensor-train object.
- `modesizes` — `ncores`-by-2 array of per-core physical row and column dimensions for the first train.

## Source

- [MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/overloads/@ttclass/sizes.m)
- [Spin Dynamics Wiki](https://spindynamics.org/wiki/index.php?title=ttclass/sizes.m)

D. Savostyanov and I. Kuprov.
