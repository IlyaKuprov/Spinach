# tests/kernel/test_polyadic_kron.m

**Source:** [https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_polyadic_kron.m](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_polyadic_kron.m)

## Purpose

Regression test for affixed polyadic tensor products and spin-space flow lifting. The test verifies that Kronecker products preserve all prefixes and suffixes, and that polyadic and sparse production generators agree.

## Behaviour

- The function creates a test result via `new_test_result` with the identifier `kernel/polyadic_kron`, the description `Affixed polyadic tensor products`, and the message `Kronecker products must preserve all prefixes and suffixes.`
- The caller's random number generator state is saved with `rng()` and restored on cleanup via `onCleanup`; the test seeds the generator with `rng(240924)`.
- Three trials build complex non-Hermitian matrices: `a` and `c` are 2-by-2 (`randn(2)+1i*randn(2)`), `b` and `d` are 3-by-3, and `ext` is a complex 2-by-3 rectangular matrix. Sparse rectangular factors are `left1` (5-by-6), `left2` (7-by-5), `right1` (6-by-4), and `right2` (4-by-8).
- A polyadic object `P=polyadic({{a,b},{c,d}})` is compared against the explicit reference `base=kron(a,b)+kron(c,d)`.
- Six object shapes are tested: `P`, `left1*P`, `left2*(left1*P)`, `P*right1`, `P*right1*right2`, and `left2*(left1*P)*right1*right2`, each against the corresponding explicit product of `base`.
- The test asserts that `objects{3}` has two prefix factors and `objects{5}` has two suffix factors.
- For each shape and both Kronecker operand orders (`kron(objects{shape},ext)` and `kron(ext,objects{shape})`), the test compares the full expansion against the explicit reference and a deferred two-column action `Q*rhs` against `ref*rhs`, where `rhs` is a complex matrix with two columns sized to the reference. Tolerances are `1e-12` absolute and relative. Labels include the trial, shape, and order indices.
- A physical one-spin Liouville system is built with `sys.magnet=14.1`, `sys.isotopes={'1H'}`, `sys.parallel={'processes',1}`, `sys.enable={'polyadic'}`, `sys.disable={'hygiene'}`, and `sys.output='hush'`; the interaction specifies `inter.zeeman.scalar={0}`. The system uses `create` and then `basis` with `bas.formalism='sphten-liouv'` and `bas.approximation='none'`. A dense comparison system is derived by clearing `sys.enable`.
- Flow parameters are `parameters.npts=[10 10]`, `parameters.dims=[0.01 0.01]`, `parameters.deriv={'period',3}`, and `parameters.v=0`.
- Three velocity cases are tested: a scalar `1e-3`, a voxel vector `1e-3*ones(10)`, and a shear `ones(10,1)*(1:10)*1e-3`.
- For each flow case, the reference is `full(v2fplanck(dense_system,parameters))` and the polyadic result is `v2fplanck(spin_system,parameters)`.
- The scalar case stores the reference as `uniform_ref`; the voxel case checks that its reference matches `uniform_ref` within `1e-12` (label `scalar versus voxel flow`).
- For non-scalar velocities, `hydrodynamics(spin_system,parameters)` is called and the test checks that `flow_x*parameters.u(:)` equals `zeros(100,1)` within `1e-12` (label `flow divergence`), confirming the transverse shear is divergence-free.
- Each flow case checks that `size(Q)` is `[400 400]` (label with ` dimensions`), that `full(Q)` matches the reference (label with ` full`), and that `Q*rhs` matches `ref*rhs` for a complex two-column right-hand side (label with ` action`), all with tolerances `1e-12` except the dimension check which uses `0`.

## Inputs and outputs

- **result** - regression results against explicit matrix references, returned by the function.
- The function takes no inputs.

## References

- [Spinach MATLAB source on GitHub](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_polyadic_kron.m)
- [Spinach documentation](https://spindynamics.org/)
