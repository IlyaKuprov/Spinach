# tests/kernel/test_giant_cache.m

Source: [https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_giant_cache.m](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_giant_cache.m)

## Purpose

Regression test that verifies the giant-spin Hamiltonian cache (`sys.enable={'ham_cache'}`) preserves giant coefficients and retention. The test checks that cached assembly returns results identical to uncached production assembly for both cache insertion orders, that changed giant coefficients are distinguished by the cache, and that the crystal decomposition receives distinct full and Zeeman Hamiltonians.

## Behaviour

- Initialises a regression record via `new_test_result` with the identifier `'kernel/giant_cache'`, the name `'Giant-spin cache identity'`, and the description `'cached assembly must preserve giant coefficients and retention.'`.
- Specifies physical Hermitian giant terms in rotated tensor frames: physical non-axial rank-two and rank-four terms have nonzero tensor and sample frames. Uncached production assembly supplies the references.
- System settings include `sys.magnet=0.73`, `sys.output='hush'`, `sys.disable={'hygiene'}`, `sys.enable={'ham_cache'}`, `sys.parallel={'processes',1}`, and `sys.parprops={}`.
- Giant interaction coefficients are `inter.giant.coeff={{[0 0 0],[2e5 3e5 7e5 -3e5 2e5],zeros(1,7),[4e4 0 2e4 0 8e4 0 2e4 0 4e4]}}` with Euler angles `inter.giant.euler={{[0 0 0],[0.31 0.53 0.17],[0 0 0],[0.19 0.41 0.37]}}`.
- Basis settings: `bas.formalism='zeeman-hilb'`, `bas.approximation='none'`. Orientation angles are `[0.29 0.47 0.61]` and isotopes are `{'E5','E6'}`.
- For each of two orderings (`ordering=1:2`), selects `sys.isotopes=isotopes(ordering)` and builds independently retained full and Zeeman-only systems via `basis(create(sys,inter),bas)` followed by `assume(spin_system,'labframe')` and `assume(spin_system,'labframe','zeeman')`.
- Obtains uncached references by clearing `sys.enable` on copies of both systems, calling `hamiltonian`, and forming `H+orientation(Q,angles)`.
- Asserts the full reference is Hermitian (`norm(references{1}-references{1}','fro')<1e-7`), has significant imaginary content (`norm(imag(references{1}),'fro')>1e5`), differs from the Zeeman-only reference (`norm(references{1}-references{2},'fro')>1e6`), and that the two references do not commute (`norm(references{1}*references{2}-references{2}*references{1},'fro')>1e12`).
- Compares both cache insertion orders and repeat cache hits by looping `n=[ordering 3-ordering ordering]`, checking `H+orientation(Q,angles)` against the corresponding reference with tolerances `1e-7` and `1e-12` under the label `'order %d retention %d'`, and checking the one-output invariant `hamiltonian(systems{n})` against the stored invariant under the label `'one-output invariant'`.
- Changes only the giant coefficients (multiplying `changed.inter.giant.coeff{1}{2}` by 1.3) with the cache still populated, verifies the changed uncached reference differs from the original by more than `1e5` in Frobenius norm, and checks the cached result against the changed reference under the label `'changed coefficients'`.
- Confirms ignored coefficients do not affect the Zeeman Hamiltonian by applying `assume(changed,'labframe','zeeman')` and checking against `references{2}` under the label `'ignored coefficients'`.
- Exercises the production crystal full and Zeeman split with `parameters.spins=sys.isotopes`, `parameters.offset=0`, `parameters.orientation=angles`, `parameters.needs={'zeeman_op'}`, calling `crystal(spin_system,@cache_generators,parameters,'labframe')` for both cached and uncached systems, and comparing under the label `'crystal decomposition'`.
- Preserves ordinary systems in both existing cache-key branches: for `modes_present=0:1`, builds a `1H` system with `inter.zeeman.scalar={2}` (and, when modes are present, adds `C3` with `inter.zeeman.scalar={2,[]}`, `inter.modes.frqs={[],1000}`, `inter.modes.anharms={[],23}`), and checks cached `hamiltonian` output against the uncached reference under the label `'ordinary cache branch'` for two repetitions.
- The local function `cache_generators` validates its arguments through `grumble` and returns `answer=[H parameters.hzeeman]` without altering the Hamiltonians.
- The local function `grumble` enforces that `spin_system` is a structure, `parameters` is a structure containing the field `hzeeman`, `H` is a numeric square matrix, and that `H`, `R`, `K`, and `parameters.hzeeman` are all numeric with matching dimensions.

## Inputs and outputs

Syntax:

```matlab
result=test_giant_cache()
```

- `result` — regression checks for retention, coefficients, and crystal, returned as the test record populated by `new_test_result` and `test_close`.
- Takes no input arguments.

## References

- Spinach GitHub source: [tests/kernel/test_giant_cache.m](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_giant_cache.m)
