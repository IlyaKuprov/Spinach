# tests/kernel/test_transform_tensor_suite.m

**Source:** [https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_transform_tensor_suite.m](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_transform_tensor_suite.m)

## Purpose

Regression test suite for the tensor transform helpers in `kernel/transform_tensor_suite`. The suite verifies that tensor transforms preserve their algebraic definitions and round-trips, covering interaction tensor parametrisations, spherical tensor round-trips, quadrupolar conversions, axial symmetrisation, and simple Hamiltonian decomposition.

## Behaviour

- Announces the test target with `fprintf('TESTING: Tensor transform helpers\n')` and initialises a test result object via `new_test_result('kernel/transform_tensor_suite', ...)`, with the failure message `'tensor transforms must preserve their algebraic definitions and round-trips.'`.
- **Haeberlen parametrisation (`anas2mat`):** with `iso=10`, `aniso=6`, `asym=0.25` and zero Euler angles, checks that the returned matrix equals `diag([iso-red_aniso*(1+asym)/2, iso-red_aniso*(1-asym)/2, iso+red_aniso])` where `red_aniso=2*aniso/3`, i.e. `diag([7, 9, 14])`, to tolerances `1e-15`.
- **Axiality/rhombicity (`axrh2mat`, `mat2axrh`):** with `iso=4`, `ax=6`, `rh=2`, checks `axrh2mat` returns `diag([2 4 6])`; `mat2axrh` on that reference returns isotropic part 4, axiality 6 (`2*zz-(xx+yy)`), rhombicity 2 (`yy-xx`), and Mehring-ordered eigenvalues `[2;4;6]` (`xx<=yy<=zz`).
- **Herzfeld-Berger span/skew (`spsk2mat`):** with `iso=5`, `span=6`, `skew=0.5`, checks the principal values `diag([1.5 6 7.5])`; span is `zz-xx` and skew fixes the middle eigenvalue displacement from isotropy.
- **Zero-field splitting (`zfs2mat`):** with `D=9`, `E=2`, checks principal values `diag([-1 -5 6])` (`[-D/3+E, -D/3-E, 2D/3]`) and tracelessness of the tensor.
- **IAS decomposition (`mat2ias`, `ias2mat`):** for `C=[1 2 -3;4 -5 6;7 8 9]`, checks the isotropic scalar equals `trace(C)/3`, the symmetric residual `A` is traceless and symmetric (`A` equals `A.'`), and that `ias2mat(a,d,A)` reconstructs `C` exactly.
- **Spherical tensor round-trip (`mat2sphten`, `sphten2mat`):** verifies that the nine spherical tensor coefficients form a complete Cartesian tensor basis by round-tripping `C`.
- **Quadratic form to spherical harmonics (`qform2sph`):** for `3*eye(3)`, checks rank-zero coefficient `6*sqrt(pi)` (`2*sqrt(pi)*a` for `Y00`), zero rank-one vector component, and zero rank-two anisotropy.
- **Stevens operators (`stev2sph`):** for `stev2sph(1,[0;1;0])`, checks the rank-one axial component equals `[0;1;0]`.
- **Traceless symmetric matrix parameters (`tsm2param`):** for `T=diag([-2 -1 3])`, checks axiality 9 (`2*3-(-2-1)`), rhombicity 1 (`yy-xx`), and that `euler2dcm(angles)*T*euler2dcm(angles).'` reconstructs `T` to optimiser tolerance (`1e-7`).
- **Quadrupolar tensor (`eeqq2nqi`):** with `Cq=1200`, `eta=0.2`, `spin_q=1`, checks principal values `diag([-240 -360 600])` — a traceless quadrupolar tensor defined by `Cq` and `eta`.
- **CASTEP EFG scaling (`castep2nqi`):** for `V=diag([-1 -2 3])`, `quad_moment=0.2`, checks the output equals `scale*V` where `scale=9.717362e+21*(quad_moment*1e-28)*1.60217657e-19/(6.62606957e-34*2*spin_q*(2*spin_q-1))`, i.e. CASTEP atomic-unit EFG tensors scaled by `e*q*Q/h/[2I(2I-1)]`.
- **WebLab cone convention (`weblab2nqi`):** with `alpha=0.11`, `theta=0.42`, `phi=0.73`, checks the first two-site tensor equals `eeqq2nqi(Cq,eta,spin_q,[-phi/2 theta alpha])` and the second equals `eeqq2nqi(Cq,eta,spin_q,[+phi/2 theta alpha])`.
- **Axial symmetrisation (`axis_tsymm`):** for `T=diag([1 3 5])` averaged around `[0;0;1]`, checks the result is `diag([2 2 5])` — full rotation around z averages xx and yy while preserving zz.
- **Hamiltonian decomposition (`ham2nqi`):** for a spin-half Hamiltonian built from Pauli operators with `omega=[11 -7 5]`, checks the recovered Zeeman vector equals `omega` and the quadrupole tensor is `zeros(3)` (spin one half has no independent quadrupolar tensor term).

## Inputs and outputs

**Syntax**

```matlab
result = test_transform_tensor_suite()
```

**Outputs**

- `result` — regression test result object with explanatory messages, accumulated through repeated calls to `test_close`.

The function takes no inputs.

## References

- Spinach source: [tests/kernel/test_transform_tensor_suite.m](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_transform_tensor_suite.m)
- Spinach GitHub repository: [https://github.com/IlyaKuprov/Spinach](https://github.com/IlyaKuprov/Spinach)
