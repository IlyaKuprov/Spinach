# tests/kernel/test_operator_expansion_suite.m

## Purpose

Regression test suite for Spinach operator expansion and conversion helpers, located at `tests/kernel/test_operator_expansion_suite.m`. It verifies that irreducible spherical tensor (IST) and bosonic monomial (BM) expansion helpers reconstruct explicit matrices, that Hilbert-to-Liouville vectorisation identities hold, and that operator-sized allocation helpers return the correct formalism dimensions.

## What the suite checks

- **Vectorisation:** for a non-symmetric two-level operator, left and right Hilbert-to-Liouville multiplication correspond to `I⊗H` and `H^T⊗I`, respectively; state-vector conversion retains column-stacking order. These identities are compared at absolute and relative tolerances of `1e-15`.
- **Operator expansions:** coefficients in irreducible spin tensors and finite bosonic monomials reconstruct the original operators, including three-level projectors, a four-level central transition and a normal-ordered boson product. The first spin energy level is counted from the bottom, whereas the first bosonic level is counted from the top. The reconstruction tolerances are `1e-12` for the IST examples and `1e-11` for the BM examples.
- **Formalism and coupling:** identity and empty-sparse allocations have dimensions `2×2` in the two-level Hilbert formalism, `4×4` in Liouville space, and `basis_dim×basis_dim` in spin-adapted Liouville space (`basis_dim=size(spin_system.bas.basis,1)`). The two-spin rank-two, zero-projection tensor is compared with `sqrt(2/3)*(Lz1*Lz2-(L+1*L-2+L-1*L+2)/4)` at `1e-15`.
- **Dissipative scaling:** for a two-level diagonal jump with requested rate `2.75`, the normalised expectation of its Lindbladian generator is checked against `-2.75` at `1e-12`.

## Inputs and outputs

**Syntax**

```matlab
result=test_operator_expansion_suite()
```

**Outputs**

- `result` — regression test result object with explanatory messages, accumulated through `test_close` and `test_true` calls.

The function takes no inputs.

## References

- Source: [tests/kernel/test_operator_expansion_suite.m](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_operator_expansion_suite.m)
