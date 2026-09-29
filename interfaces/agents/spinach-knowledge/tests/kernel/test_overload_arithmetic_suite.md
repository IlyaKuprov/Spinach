# tests/kernel/test_overload_arithmetic_suite.m

**Source:** [https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_overload_arithmetic_suite.m](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_overload_arithmetic_suite.m)

## Purpose

Regression test for the "cheap" overload arithmetic in Spinach covering cell arrays, structs, RCV sparse storage, and polyadic classes. The test verifies that these overloaded operations match explicit Matlab matrix arithmetic on small examples.

## Behavior

The function announces the test target with `fprintf('TESTING: Cheap overload arithmetic\n')` and initializes a test result object via `new_test_result` with suite name `kernel/overload_arithmetic_suite`, description `Cheap overload arithmetic`, and the message `cell, struct, RCV, and polyadic overloads must match explicit Matlab arithmetic on small examples.`

Each check is performed with `test_close`, which appends pass/fail results and explanatory messages to the returned result object.

### Cell overloads

Using operands `A=[1 2;3 4]`, `B=[0 5;-1 2]`, `cell_a={A,B}`, and `cell_b={eye(2),ones(2)}`, the test checks:

- Cell plus: `cell_a+cell_b` adds corresponding cells element by element, compared against `A+eye(2)` and `B+ones(2)`.
- Cell minus: `cell_a-cell_b` subtracts corresponding cells, compared against `A-eye(2)`.
- Scalar plus on the right (`cell_a+2`) and left (`2+cell_a`), applying the scalar to every cell.
- Elementwise times with a scalar array: `cell_a.*[2 3]` and `[2 3].*cell_a` multiply each cell by the matching scalar (`2*A`, `3*B`).
- Matrix multiplication with `R=diag([2 3])`: `cell_a*R` and `R*cell_a` multiply every cell from the right or left respectively.
- `totsum` on a cell of sparse matrices `{sparse([1 0;0 2]),sparse([0 3;4 0])}` yields `sparse([1 3;4 2])`, preserving sparse arithmetic.
- `inflate(cell_a)` leaves each numeric entry unchanged.
- `complex(cell_a)` applies `complex` to every cell.

All cell checks use tolerances `1e-15` for both absolute and relative error.

### Struct overloads

With structs `s1.alpha=[1;2]`, `s1.beta.gamma=[3 4]`, `s2.alpha=[5;6]`, `s2.beta.gamma=[7 8]`:

- `s1+s2` recursively adds matching numeric fields: `s3.alpha` equals `[6;8]` and `s3.beta.gamma` equals `[10 12]`.
- `2*s1` applies numeric left multiplication recursively: `s4.alpha` equals `[2;4]` and `s4.beta.gamma` equals `[6 8]`.

Tolerances are `1e-15`.

### RCV sparse storage

Using `S=sparse([1 0 2;0 3 0])`, `T=sparse([0 4 0;5 0 6])`, `R1=rcv(S)`, `R2=rcv(T)`, and `C=sparse([1+1i 0;0 2-3i])`, the test checks:

- `full(R1)` reconstructs the original sparse matrix.
- `R1+R2` matches `S+T` (RCV plus concatenates entries; sparse conversion sums duplicates).
- `R1+T` matches `S+T` (RCV plus Matlab sparse).
- `R1-R2` matches `S-T`.
- `2*R1` matches `2*S` (scalar-matrix multiplication scales stored values).
- `3.*R1` matches `3*S` (scalar-array multiplication scales stored values).
- `R1./2` matches `S./2` (right division by a scalar scales stored values).
- `R1.'` matches `S.'` (transpose swaps stored row and column indices).
- `rcv(C)'` matches `C'` (conjugate transpose swaps indices and conjugates values).
- `R1*rcv(S.')` matches `S*S.'` (matrix multiplication).
- `[R1 R2]` matches `[S T]` (horizontal concatenation shifts right-hand column indices).
- `[R1;R2]` matches `[S;T]` (vertical concatenation shifts lower row indices).
- `size(R1)` reports `[2 3]` as stored matrix dimensions.

All RCV checks use tolerances `1e-15`.

### Polyadic arithmetic

Using cores `P1=[1 2;0 3]`, `P2=[2 -1;4 1]`, `P3=[0 1;5 2]`, `P4=[3 0;-2 1]`, the polyadic object `P=polyadic({{P1,P2},{P3,P4}})` is compared against the reference `P_ref=kron(P1,P2)+kron(P3,P4)`:

- `full(P)` opens and sums all stored Kronecker products (tolerance `1e-14`).
- `size(P)` equals `[4 4]`, the product of core dimensions (tolerance `1e-15`).
- `nnz(P)` equals `nnz(P1)+nnz(P2)+nnz(P3)+nnz(P4)`, counting non-zero entries in all cores (tolerance `1e-15`).
- `full(P+speye(4))` matches `P_ref+eye(4)`: adding a matrix buffers a new additive term (tolerance `1e-14`).
- `full(2*P)` matches `2*P_ref`: left scalar multiplication scales the buffered cores (tolerance `1e-14`).
- `full(P*3)` matches `3*P_ref`: right scalar multiplication scales the buffered cores (tolerance `1e-14`).
- With `v=(1:4)'`, `P*v` matches `P_ref*v`: multiplication by a dense vector matches the opened matrix (tolerance `1e-14`).
- `full(P.')` matches `P_ref.'`: transpose transposes every core and swaps prefixes with suffixes (tolerance `1e-14`).
- `full(P')` matches `P_ref'`: conjugate transpose conjugates and transposes every core (tolerance `1e-14`).
- `full(kron(polyadic({{P1}}),P2))` matches `kron(P1,P2)`: polyadic kron appends matrix cores without opening the product (tolerance `1e-14`).

## Inputs and outputs

```matlab
result = test_overload_arithmetic_suite()
```

- **Output:** `result` — regression test result object with explanatory messages, produced by `new_test_result` and accumulated through repeated `test_close` calls.
- **Input:** none.

## References

- [Spinach GitHub repository — test_overload_arithmetic_suite.m](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_overload_arithmetic_suite.m)
- [Spinach project website](http://spindynamics.org/)
