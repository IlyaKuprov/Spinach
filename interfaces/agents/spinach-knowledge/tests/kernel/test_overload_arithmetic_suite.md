# tests/kernel/test_overload_arithmetic_suite.m

**Source:** [https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_overload_arithmetic_suite.m](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_overload_arithmetic_suite.m)

## Purpose

Regression test for the "cheap" overload arithmetic in Spinach covering cell arrays, structs, RCV sparse storage, and polyadic classes. The test verifies that these overloaded operations match explicit Matlab matrix arithmetic on small examples.

## Arithmetic invariants

- **Cell arrays:** elementwise addition and subtraction, scalar addition from either side, matched elementwise scaling, and left/right matrix multiplication act on each numeric cell as ordinary matrix arithmetic would. Summing sparse cells retains the expected sparse matrix; `inflate` leaves their values unchanged and `complex` converts each cell. These comparisons use `1e-15` absolute and relative tolerance.
- **Structs:** addition and matrix multiplication distribute recursively over numeric fields, including a nested field, and agree with explicit field arithmetic at `1e-15`.
- **RCV sparse storage:** `full()` reconstructs the sparse reference, including sums of duplicate stored coordinates. RCV–RCV and RCV–sparse addition, subtraction, scalar scaling/division, ordinary and conjugate transpose, matrix multiplication, horizontal/vertical concatenation, and size all agree with standard sparse matrices at `1e-15`.
- **Polyadic sums of Kronecker products:** the stored polyadic object opens to `P_ref=kron(P1,P2)+kron(P3,P4)`. Its size is `[4 4]`, while `nnz` counts nonzero entries in the stored cores. Matrix addition, scalar and vector multiplication, ordinary/conjugate transpose, and appended Kronecker factors agree with their explicit counterparts; the matrix comparisons use `1e-14` tolerance, with size and core-count checks at `1e-15`.

## Inputs and outputs

```matlab
result=test_overload_arithmetic_suite()
```

- **Output:** `result` — regression test result summarising the algebraic identity checks and any failed assertions.
- **Input:** none.

## References

- [Spinach GitHub repository — test_overload_arithmetic_suite.m](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_overload_arithmetic_suite.m)
- [Spinach project website](http://spindynamics.org/)
