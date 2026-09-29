# tests/kernel/test_dynamic_frame_frontends.m

Source: https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_dynamic_frame_frontends.m

## Purpose

Regression test for the frame-transformation and averaging front-end kernels. It exercises `carrier()`, `frqoffset()`, `rotframe()`, `average()`, and `orientation()` on compact systems and analytically controlled limiting cases, checking that the frame-transformation helpers preserve exact algebraic limits.

## Behaviour

The test announces its target with `fprintf('TESTING: Dynamic frame and averaging front ends\n')` and initialises a result object via `new_test_result('kernel/dynamic_frame_frontends', ...)`. It then runs five local subtests:

- **Carrier algebra** (`local_test_carrier`): builds a heteronuclear Liouville-space system (`sphten-liouv` formalism, `1H`/`13C`, 14.1 T magnet, scalar Zeeman 1.0 and 2.0, scalar coupling 10.0 between spins 1 and 2) and evaluates `carrier()` in `left`, `right`, `comm`, and `acomm` modes for the proton channel. It checks that the `comm` mode equals `left` minus `right` side products and that `acomm` equals `left` plus `right`, both with tolerances `1e-8` (relative) and `1e-14` (absolute). It also checks that `carrier(spin_system,'all','comm')` equals the sum of the `1H` and `13C` comm-mode carriers, and that the proton `comm` carrier equals `basefrqs(1)` times `operator(spin_system,{'Lz'},{1},'comm')`.
- **Frequency offsets** (`local_test_frqoffset`): starting from a zero Hamiltonian, applies `frqoffset()` with `parameters.spins={'1H','13C'}` and `parameters.offset=[25.0,-40.0]`, comparing against `2*pi*25.0*operator(spin_system,'Lz','1H') - 2*pi*40.0*operator(spin_system,'Lz','13C')` with tolerances `1e-12`/`1e-12`. A second call with duplicated proton channels (`{'1H','13C','1H'}` and offsets `[25.0,-40.0,25.0]`) must produce the same Hamiltonian, verifying that duplicate channel names with equal offsets are merged.
- **Rotating frame** (`local_test_rotframe`): builds a one-spin Hilbert-space laboratory-frame system (`zeeman-hilb` formalism, `1H`, 14.1 T, zero scalar Zeeman, then `assume(spin_system,'labframe')`). Forms `H0=carrier(spin_system,'1H')` and a symmetrised transverse perturbation `H1=1e3*operator(spin_system,'Lx','1H')`. Checks that `rotframe(spin_system,H0,H0+H1,'1H',0)` returns `H1` with tolerances `1e-10`/`1e-14`, i.e. the zeroth-order transformation removes only the carrier.
- **Averaging theories** (`local_test_average`): defines an unmodulated two-level decomposition with `Hp` and `Hm` zero and `H0=sparse([0 1;1 0])`, at `omega=2*pi*1000`. Loops over the theories `{'ah_first_order','ah_second_order','ah_third_order','kb_first_order','kb_second_order','kb_third_order','matrix_log'}` and verifies each returns `H0` unchanged (tolerances `1e-10`/`1e-12`).
- **Orientation contraction** (`local_test_orientation`): builds a hand-made sparse rotational basis `Q` with rank-one (3-by-3) and rank-two (5-by-5) blocks of 2-by-2 sparse matrices, populating `Q{1}{1,3}`, `Q{1}{3,1}`, and `Q{2}{3,3}`. Contracts at Euler angles `[0.2 0.3 -0.4]` via `orientation(Q,angles)` and compares against the Wigner-matrix combination `D1(1,3)*Q{1}{1,3}+D1(3,1)*Q{1}{3,1}+D2(3,3)*Q{2}{3,3}` (symmetrised), with tolerances `1e-14`/`1e-14`.

Each check is recorded through `test_close(result,name,observed,reference,rtol,atol,message)`, which appends regression results with explanatory messages.

## Inputs and outputs

Syntax:

```matlab
result=test_dynamic_frame_frontends()
```

- **Output**: `result` — regression test result object with explanatory messages, accumulated across all subtests.
- **Input**: none.

## References

- Source file: https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_dynamic_frame_frontends.m
