# tests/kernel/test_operator_expansion_suite.m

## Purpose

Regression test suite for Spinach operator expansion and conversion helpers, located at `tests/kernel/test_operator_expansion_suite.m`. It verifies that irreducible spherical tensor (IST) and bosonic monomial (BM) expansion helpers reconstruct explicit matrices, that Hilbert-to-Liouville vectorisation identities hold, and that operator-sized allocation helpers return the correct formalism dimensions.

## Behavior

The function announces the test target with `fprintf`, initializes a result object via `new_test_result` keyed to `kernel/operator_expansion_suite` with the message that operator expansion coefficients must reconstruct the source matrices, and then runs a sequence of checks:

- **Hilbert-to-Liouville vectorisation** on the non-diagonal matrix `H=[1 2;3 4]` with `unit=speye(2)`: `hilb2liouv(H,'left')` against `kron(unit,H)`; `hilb2liouv(H,'right')` against `kron(transpose(H),unit)`; `hilb2liouv(H,'statevec')` against `H(:)`. Tolerances are `1e-15` for both absolute and relative closeness.
- **IST expansion** of the generic spin-one matrix `A=[1 2 0;0 -1 3;4 0 2]` via `oper2ist`, reconstructed with the local `ist_reconstruct` helper against `A` at tolerances `1e-12`.
- **Energy-level counting conventions** via `enlev2ist(3,1,'S')` reconstructed against `P=zeros(3); P(3,3)=1` (spin levels counted from the bottom upward), and `enlev2ist(3,1,'B')` against `P=zeros(3); P(1,1)=1` (boson levels counted from the top downward), both at `1e-12` tolerances.
- **Central-transition wrapper** `ct2ist(4,'+')` reconstructed against `centrans(4,'+')` at `1e-12` tolerances.
- **Boson-product wrapper** `bos2ist('CAN',3)` reconstructed against `W.c*W.a*W.n`, where `W=weyl(3)`, at `1e-12` tolerances.
- **Bosonic monomial expansion** of the same generic matrix `A` via `oper2bm` at `1e-11` tolerances, and `enlev2bm(3,2)` reconstructed against `P=zeros(3); P(2,2)=1` at `1e-11` tolerances.
- **Formalism dimensions** using a one-spin system (`sys.magnet=0`, `sys.isotopes={'1H'}`, `inter.zeeman.scalar={0}`, `bas.approximation='none'`) built with `test_spin_system`:
  - `zeeman-hilb`: `unit_oper` against `speye(2)` at `1e-15`; `mprealloc(spin_system,2)` must be an empty (`nnz(A)==0`) sparse `2x2` matrix.
  - `zeeman-liouv`: `unit_oper` against `speye(4)` at `1e-15`; `mprealloc(spin_system,2)` must be an empty sparse `4x4` matrix.
  - `sphten-liouv`: `unit_oper` against `speye(basis_dim)` where `basis_dim=size(spin_system.bas.basis,1)`; `mprealloc(spin_system,2)` must be an empty sparse `basis_dim x basis_dim` matrix.
- **Two-spin irreducible tensor formula** in Hilbert space with `sys.isotopes={'1H','1H'}`, `inter.zeeman.scalar={0,0}`, `inter.coupling.scalar{1,2}=0`, `inter.coupling.scalar{2,2}=0`, `bas.formalism='zeeman-hilb'`: `twospinist(spin_system,1,2,[2 0],'comm')` is compared against the reference `sqrt(2/3)*(operator(spin_system,{'Lz','Lz'},{1,2})-(1/4)*(operator(spin_system,{'L+','L-'},{1,2})+operator(spin_system,{'L-','L+'},{1,2})))` at `1e-15` tolerances.
- **Lindbladian rate calibration** on a diagonal jump process with `rho=[1;1]` and `rate=2.75`: `R=lindbladian(diag([1 0]),diag([0 1]),rho,rate)` and the observable `obs=real((rho'*R*rho)/(rho'*rho))` is compared against `-rate` at `1e-12` tolerances.\n
Local helper functions:

- `ist_reconstruct(mult,states,coeffs)` builds the complete IST basis with `irr_sph_ten(mult)` and recombines zero-based Spinach IST indices into a matrix as `A=A+coeffs(n)*T{states(n)+1}` over all states.
- `bm_reconstruct(nlevels,states,coeffs)` builds the complete bosonic monomial basis with `boson_mono(nlevels)` and recombines zero-based Spinach BM indices analogously.

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
