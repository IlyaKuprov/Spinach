# tests/kernel/test_indexing_roundtrip_suite.m

- Signature: `result=test_indexing_roundtrip_suite()`

## Purpose

Regression test for angular-momentum and matrix indexing helpers. It checks that linear and structured indices are mutually consistent for spherical tensors, Wigner functions, and serpentine matrix indexing.

## Checks

1. **Spin-state round trip:** For zero-based linear indices `I=0:24`, `lin2lm(I)` produces `L,M`, and `lm2lin(L,M)` must equal `I` exactly.
2. **Spin-state ordering:** For the first nine indices, `L(1:9)` must be `[0 1 1 1 2 2 2 2 2]` and `M(1:9)` must be `[0 1 0 -1 2 1 0 -1 -2]`: ranks increase, with projections decreasing within each rank.
3. **Wigner-function round trip:** For one-based linear indices `J=1:35`, `lin2lmn(J)` produces `Lw,Mw,Nw`, and `lmn2lin(Lw,Mw,Nw)` must equal `J` exactly.
4. **Wigner rank-one ordering:** `[Lw(2:10); Mw(2:10); Nw(2:10)]` must equal `[ones(1,9); 1 1 1 0 0 0 -1 -1 -1; 1 0 -1 1 0 -1 1 0 -1]`: ranks increase, then `M` decreases, then `N` decreases within each `M`.
5. **Serpentine tables:** `S1=serpentine(4,1)` must equal `[1 3 6 10; 2 5 9 13; 4 8 12 15; 7 11 14 16]`, the base-one anti-diagonal triangular scan order. `S0=serpentine(4,0)` must equal that table minus one.
6. **Serpentine coordinate round trips:** With `N=4`, `lin2kq(N,idx1,1)` followed by `kq2lin(N,K1,Q1,1)` must reproduce `idx1=1:N^2`; `lin2kq(N,idx0,0)` followed by `kq2lin(N,K0,Q0,0)` must reproduce `idx0=0:(N^2-1)`. Both must be exact.

All comparisons use `test_close` with zero tolerances. The test announces `TESTING: Indexing conversion functions` and initializes its result with `new_test_result('kernel/indexing_roundtrip_suite','Indexing conversion functions','indexing helpers must be exact inverses on valid integer domains.')`.

## Outputs

- `result` — regression test result with explanatory messages.