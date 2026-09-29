# tests/kernel/test_indexing_roundtrip_suite.m

## Purpose

Regression test for Spinach's angular-momentum and matrix indexing helper functions. It verifies that linear and structured index representations are mutually consistent for spherical tensors, Wigner functions, and matrix serpentine indexing.

## Behaviour

The function announces the test target with `fprintf('TESTING: Indexing conversion functions\n')` and initialises a test result object via `new_test_result` with the identifier `kernel/indexing_roundtrip_suite`, the description `Indexing conversion functions`, and the criterion `indexing helpers must be exact inverses on valid integer domains.`

It then performs the following checks, each through `test_close` with zero tolerances (exact comparison):

- **`lin2lm`/`lm2lin` round-trip**: for `I=0:24`, `lm2lin(L,M)` must return `I` exactly.
- **`lin2lm` first ranks**: `L(1:9)` must equal `[0 1 1 1 2 2 2 2 2]`, i.e. linear LM state indexing lists rank 0, then the three rank-1 states, then rank 2.
- **`lin2lm` first projections**: `M(1:9)` must equal `[0 1 0 -1 2 1 0 -1 -2]`, i.e. within each rank, projections are listed in decreasing `M` order.
- **`lin2lmn`/`lmn2lin` round-trip**: for `J=1:35`, `lmn2lin(Lw,Mw,Nw)` must return `J` exactly.
- **`lin2lmn` rank-one block**: `[Lw(2:10); Mw(2:10); Nw(2:10)]` must equal `[ones(1,9); 1 1 1 0 0 0 -1 -1 -1; 1 0 -1 1 0 -1 1 0 -1]`, i.e. rank-one Wigner functions are ordered by decreasing `M` and, within each `M`, decreasing `N`.
- **`serpentine` base one**: `serpentine(4,1)` must equal `[1 3 6 10; 2 5 9 13; 4 8 12 15; 7 11 14 16]` (anti-diagonal triangular ordering).
- **`serpentine` base zero**: `serpentine(4,0)` must equal the base-one reference minus one.
- **`lin2kq`/`kq2lin` base one**: for `N=4` and `idx1=1:N^2`, `kq2lin(N,K1,Q1,1)` must return `idx1` exactly.
- **`lin2kq`/`kq2lin` base zero**: for `idx0=0:(N^2-1)`, `kq2lin(N,K0,Q0,0)` must return `idx0` exactly.

Linear spin-state indexing is zero-based and ordered by increasing `L`; Wigner D-function indexing is one-based and ordered by increasing `L`, then `M`, then `N`.

## Inputs and outputs

- **Output**: `result` — regression test result with explanatory messages, accumulated by successive `test_close` calls.
- **Input**: none.

## References

- Source: [tests/kernel/test_indexing_roundtrip_suite.m](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_indexing_roundtrip_suite.m)
