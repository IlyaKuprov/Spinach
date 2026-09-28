# tests/kernel/test_indexing_inverse_suite.m

- Signature: `result=test_indexing_inverse_suite()`

## Checks

- Compares `serpentine(4,1)` with the base-one table `[1 3 6 10; 2 5 9 13; 4 8 12 15; 7 11 14 16]` and `serpentine(4,0)` with that table minus one.
- Checks `kq2lin` against both 4×4 tables and checks that `lin2kq` recovers every row and column index, using `1:4` for base one and `0:3` for base zero.
- For every `(L,M)` with `L=0:5` and `M` descending from `L` to `-L`, checks that `lm2lin` produces consecutive zero-based indices in that order and that `lin2lm` recovers both coordinates.
- For every `(L,M,N)` with `L=0:3` and both `M` and `N` descending from `L` to `-L`, checks that `lmn2lin` produces consecutive one-based indices in that order and that `lin2lmn` recovers all three coordinates.

## Output

- `result` — regression test result with explanatory messages.