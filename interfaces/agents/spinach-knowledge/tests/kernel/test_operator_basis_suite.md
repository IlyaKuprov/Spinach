# tests/kernel/test_operator_basis_suite.m

## Purpose

Regression test suite for operator-basis construction and expansion helpers in Spinach. It verifies that irreducible spherical tensors, Stevens operators, Weyl boson operators, bosonic monomials, single-transition and central-transition operators, and the various expansion/reconstruction helpers satisfy their defining algebra and reproduce expected operators.

Source: [tests/kernel/test_operator_basis_suite.m](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_operator_basis_suite.m)

## Behaviour

- Announces the test target with `fprintf('TESTING: Operator-basis construction functions\n')` and initialises a test result object via `new_test_result('kernel/operator_basis_suite', ...)`.
- **Irreducible spherical tensors**: for multiplicity 3, checks `[Lz, T(k,m)] = m*T(k,m)` for rank-2 tensors with projections `2:-1:-2` (tolerance `1e-13`), and that `irr_sph_ten(mult)` returns `mult^2` operators.
- **Stevens operators**: verifies `stevens(3,1,0)` equals the angular momentum `Lz` operator (tolerance `1e-14`).
- **Weyl boson operators**: for `weyl(4)`, checks `c*a = n`, `[n,c] = c`, and `[n,a] = -a` (tolerance `1e-14`).
- **Bosonic monomials**: checks `boson_mono(3)` first element is the identity (`speye(3)`) and second element is the creation operator `weyl(3).c`; checks `boson_ortho(3)` elements have zero mutual Hilbert-Schmidt overlap (tolerance `1e-12`).
- **Single-transition operators**: for `sin_tran(3)`, checks each of the 9 matrices equals a single-entry sparse matrix at the serpentine index location given by `lin2kq(3,n,1)`.
- **Central-transition operators**: checks `centrans(4,'z')` equals `sparse([2 3],[2 3],[0.5 -0.5],4,4)` and `centrans(4,'+')` equals `sparse(2,3,1,4,4)` (tolerance `1e-14`).
- **Expansion helpers**: verifies reconstruction (tolerance `1e-13`) for:
  - `oper2ist` on a 2x2 complex matrix via `ist_reconstruct`.
  - `ct2ist(4,'z')` reconstructing `centrans(4,'z')`.
  - `enlev2ist(3,2,'S')` reconstructing the projector onto level 2.
  - `bos2ist('CA',3)` reconstructing `weyl(3).c*weyl(3).a`.
  - `oper2bm` and `enlev2bm(3,2)` reconstructing the same projector via `bm_reconstruct`.
- **Two-spin IST**: builds a two-proton spin system (`sys.magnet=0`, isotopes `{'1H','1H'}`, Zeeman-Hilbert formalism) and checks `twospinist(spin_system,1,2,[2 0],'comm')` against the Cartesian product expression `sqrt(2/3)*(Lz*Lz - (L+*L- + L-*L+)/4)` (tolerance `1e-14`).
- **Sparse preallocation**: checks `mprealloc(spin_system,2)` returns size `[4 4]` under `zeeman-hilb` formalism and `[16 16]` under `zeeman-liouv` formalism.

## Inputs and outputs

**Syntax**

```matlab
result = test_operator_basis_suite()
```

**Outputs**

- `result` — regression test result object with explanatory messages for each check.

The function takes no inputs.

## References

- [Spinach GitHub repository — test_operator_basis_suite.m](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_operator_basis_suite.m)
