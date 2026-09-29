# tests/kernel/test_rotor_assume.m

Source: [https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_rotor_assume.m](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_rotor_assume.m)

## Purpose

Regression test for explicit assumptions during rotor-stack construction. It verifies that explicit assumptions reach numerical rotating frames, that fresh, stale, and matching objects produce identical numerical-frame stacks, and that inconsistent numerical-frame requests are rejected.

## Behaviour

- The test is registered via `new_test_result` as `kernel/rotor_assume` with the description `Rotor-stack assumptions` and the message `Explicit assumptions must reach numerical rotating frames.`
- The first system is an anisotropic heteronuclear pair (`1H`, `13C`) at 9.4 T with noncommuting laboratory Hamiltonians: Zeeman eigenvalues `[-12 5 20]` and `[-30 10 45]`, Zeeman Euler angles `[0.2 0.5 0.7]` and `[0.6 0.4 0.3]`, and a scalar coupling of 150 between the two spins. The rotor axis is `[1 2 3]/sqrt(14)`, transmitter offsets are `[70 -35]`, `max_rank` is 1, and the orientation is `[0.3 0.7 0.2]`.
- Both formalisms `zeeman-hilb` and `sphten-liouv` are exercised, together with both MAS orientation conventions `rotor` and `magnet` (`masframe`), and rotating frames `{{'1H',1},{'13C',1}}`.
- For each formalism and MAS frame, the test compares rotor stacks from a matching `labframe` assumption, a fresh spin system, and a stale object carrying a prior `nmr` assumption, using a tolerance `frame_tol = 1e-14*carrier_norm*sqrt(2*parameters.max_rank+1)` where `carrier_norm` is the sum of Frobenius norms of the `1H` and `13C` carriers. Repeated calls on the matching object must reproduce the reference stack. Rotor phase grids must be identical (tolerances 0) regardless of assumption history.
- Requests for `nmr` numerical frames on spins already in the rotating frame are expected to fail with an error message containing `already in the rotating frame`.
- Empty-frame (`rframes={}`) NMR rotor stacks must be independent of prior assumptions (compared with tolerances 0).
- The laboratory control stack must be complex (`norm(imag(lab{1}),'fro')>1`) and noncommuting (`norm(lab{1}*lab{2}-lab{2}*lab{1},'fro')>1`).
- The second system is a hyperfine-coupled pair (`E`, `1H`) at 0.35 T with electron Zeeman eigenvalues `[2.0023 2.0023 2.0023]`, nuclear Zeeman eigenvalues `[-12 5 20]`, hyperfine coupling eigenvalues `[1e4 2e4 4e4]`, and hyperfine Euler angles `[0.4 0.6 0.8]`. The `masframe` is `rotor`.
- Assumption sets tested are `nmr`, `cavity`, `esr`, `deer`, `deer-zz`, and `spin-phonon`. For each, empty-frame rotor stacks from fresh and assumed objects must match (tolerances 0).
- Direct `rotframe` transformations of already-rotating spins are rejected with `already in the rotating frame`. Targets are `E` only, except for `nmr` and `cavity` where both `E` and `1H` are targeted. The same rejection is required for `rotor_stack` calls with explicit rotating frames on both fresh and stale (`labframe`-assumed) inputs.
- Nuclear numerical frames (`rframes={{'1H',1}}`) under the electron-only assumptions `esr`, `deer`, `deer-zz`, and `spin-phonon` must match the `esr` reference within `frame_tol = 1e-14*norm(carrier(spin_system,'1H'),'fro')*sqrt(numel(fresh))`, for both fresh and stale inputs.
- The mixed-frame control (`spin-phonon`, empty frames) must be complex and noncommuting, attributed to pseudosecular hyperfine terms.
- The third configuration zeroes the hyperfine coupling (`[0 0 0]`) and offsets (`[0 0]`) and tests the solid-effect components `se_dnp_h+`, `se_dnp_h-`, and `se_dnp_h0`. Empty-frame rotor stacks must contain only zero matrices. Direct `rotframe` and `rotor_stack` numerical-frame requests for both `E` and `1H` must be rejected with an error message containing `solid-effect Hamiltonian components`.
- The fourth system pairs `1H` with the quadrupolar nucleus `14N` at 9.4 T, with `14N` quadrupolar coupling eigenvalues `[-1e5 -2e5 3e5]` and Euler angles `[0.4 0.6 0.8]`. Under the `qnmr` assumption, numerical frames on `14N` (`rframes={{'14N',1}}`) must match the reference within `frame_tol = 1e-14*norm(carrier(assumed,'14N'),'fro')*sqrt(numel(fresh))` for both fresh and stale (`nmr`-assumed) inputs, while a `1H` numerical frame under `qnmr` must be rejected with `already in the rotating frame`.

## Inputs and outputs

```matlab
result=test_rotor_assume()
```

- **Output**: `result` — regression test result object with explanatory messages, accumulated through `test_close` and `test_true` checks.
- The function takes no inputs.

## References

- Source file: [tests/kernel/test_rotor_assume.m](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_rotor_assume.m) in the Spinach repository.
