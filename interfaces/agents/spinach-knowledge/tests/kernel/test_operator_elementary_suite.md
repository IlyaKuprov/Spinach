# tests/kernel/test_operator_elementary_suite.m

**Source:** [https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_operator_elementary_suite.m](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_operator_elementary_suite.m)

## Purpose

Regression test for the elementary operator generators in `kernel/operators`. The suite verifies analytic commutation relations, indexing conventions, and small explicit matrices for the low-level operator constructors.

## Behaviour

The function announces the test target with `fprintf('TESTING: Elementary operator generators\n')` and initialises a test result object via `new_test_result('kernel/operator_elementary_suite', ...)`, describing the requirement that low-level operator constructors satisfy their defining algebraic identities. Each subsequent check is performed with `test_close`, which compares a computed quantity against a reference within specified absolute and relative tolerances and appends explanatory messages to the result.

The checks performed are:

- **`pauli(3)`** — spin-one angular momentum: `[Sx,Sy] = i Sz` (tolerances `1e-14`); the raising operator `S.p` equals `S.x + 1i*S.y` (`1e-15`); the unit operator `S.u` equals `speye(3)` (`1e-15`).
- **`weyl(4)`** — finite-truncation Weyl algebra, checked away from the unavoidable edge state: `W.c*W.a` equals the number operator `W.n` (`1e-14`); `[W.n, W.c] = W.c` (`1e-14`); `[W.n, W.a] = -W.a` (`1e-14`).
- **`boson_mono(3)`** — bosonic monomial serpentine indexing: `B{1}` equals `W.u(1:3,1:3)` for the `weyl(4)` unit (`1e-15`); `B{2}` equals `weyl(3).c` (`1e-15`); `B{3}` equals `weyl(3).a` (`1e-15`).
- **`boson_ortho(3)`** — Gram–Schmidt orthogonality without imposing normalisation: the Gram matrix of Hilbert–Schmidt inner products `trace(full(B{n}'*B{k}))` is compared against its own diagonal part, i.e. off-diagonal overlaps must vanish (`1e-12`).
- **`sin_tran(4)`** — single-transition basis indexing from the documented 4×4 map: serpentine index ten `A{10}` equals `sparse(1,4,1,4,4)` (`1e-15`); serpentine index seven `A{7}` equals `sparse(4,1,1,4,4)` (`1e-15`).
- **`centrans(4, ...)`** — central-transition operators embedded in a spin-3/2 manifold: `centrans(4,'z')` equals a sparse 4×4 matrix with entries `ct_z(2,2)=0.5`, `ct_z(3,3)=-0.5` (`1e-15`); `centrans(4,'+')` equals a matrix with `ct_p(2,3)=1` (`1e-15`); `centrans(4,'-')` equals a matrix with `ct_m(3,2)=1` (`1e-15`).
- **`irr_sph_ten(3,2)`** — irreducible spherical tensor projection quantum numbers: for projections `proj = [2 1 0 -1 -2]` and `L = pauli(3)`, each tensor obeys `[L.z, T{n}] = proj(n)*T{n}` (`1e-12`).
- **`stevens(3, ...)`** — Stevens operators: `stevens(3,0,0)` equals `speye(3)` (`1e-15`); `stevens(3,2,+1)` and `stevens(3,2,-1)` are each Hermitian, checked against their own conjugate transposes (`1e-14`).

## Inputs and outputs

```matlab
result = test_operator_elementary_suite()
```

- **`result`** — regression test result object with explanatory messages, accumulated from the individual `test_close` checks.

The function takes no inputs.

## References

- [Spinach GitHub repository](https://github.com/IlyaKuprov/Spinach)
