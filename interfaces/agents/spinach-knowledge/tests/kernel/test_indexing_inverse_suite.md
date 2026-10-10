# tests/kernel/test_indexing_inverse_suite.m

**Source:** [https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_indexing_inverse_suite.m](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_indexing_inverse_suite.m)

## Purpose

Regression test suite that verifies the indexing helper functions in Spinach are exact inverses of each other on their finite integer domains. The suite checks serpentine matrix indexing (`serpentine`, `kq2lin`, `lin2kq`), spin-state L,M indexing (`lm2lin`, `lin2lm`), and Wigner-function L,M,N indexing (`lmn2lin`, `lin2lmn`) over complete low-rank domains.

## Behaviour

- Announces the test target with `TESTING: Indexing inverse helpers` and initialises a test result object named `kernel/indexing_inverse_suite` with the description "Indexing inverse helpers" and the requirement that indexing helpers must be exact inverses on their finite integer domains.
- **Serpentine base-one check:** compares `serpentine(4,1)` against the documented reference matrix `[1 3 6 10;2 5 9 13;4 8 12 15;7 11 14 16]`, confirming that base-one serpentine indexing follows anti-diagonals from lower to upper rows.
- **Serpentine base-zero check:** compares `serpentine(4,0)` against the base-one reference matrix minus one, confirming the base-zero table is the base-one table shifted by one.
- **k,q round-trip, base one:** builds `K` and `Q` grids over `1:4`, converts with `kq2lin(4,K,Q,1)`, inverts with `lin2kq(4,I,1)`, and verifies both the returned row and column indices survive the round trip; also checks that the linear index table equals the base-one serpentine reference matrix.
- **k,q round-trip, base zero:** repeats the same checks with `K` and `Q` grids over `0:3`, `kq2lin(4,K,Q,0)`, and `lin2kq(4,I,0)`, verifying row and column round trips and that the linear index table equals the base-one reference matrix minus one.
- **L,M ordering and inverse:** constructs a complete low-rank domain with `l` from 0 to 5 and, for each rank, `m` running from `l` down to `-l`; verifies that `lm2lin(L,M)` produces linear indices `0:(numel(I)-1)` (indices increase rank by rank and projection by decreasing M), and that `lin2lm(I)` recovers every L rank and M projection exactly.
- **L,M,N ordering and inverse:** constructs a complete Wigner-function domain with `l` from 0 to 3, `m` from `l` down to `-l`, and `n` from `l` down to `-l` for each (l,m) pair; verifies that `lmn2lin(L,M,N)` produces linear indices `1:numel(I)` (indices start at one and then increase rank by rank), and that `lin2lmn(I)` recovers every L rank, M projection, and N projection exactly.
- All comparisons use `test_close` with zero tolerances, so the helpers must match exactly.

## Inputs and outputs

```matlab
result = test_indexing_inverse_suite()
```

- **Outputs**
  - `result` — regression test result object with explanatory messages for each check.
- **Inputs**
  - None.

## References

- `serpentine` — serpentine matrix indexing helper.
- `kq2lin`, `lin2kq` — k,q to linear index conversion and its inverse.
- `lm2lin`, `lin2lm` — spin-state L,M to linear index conversion and its inverse.
- `lmn2lin`, `lin2lmn` — Wigner-function L,M,N to linear index conversion and its inverse.
- `new_test_result`, `test_close` — test harness utilities.
