# tests/kernel/test_hilb_pulse_prop.m

**Source:** [https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_hilb_pulse_prop.m](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_hilb_pulse_prop.m)

## Purpose

Regression test for Hilbert-space shaped-pulse propagators and density-state reuse, exercising `shaped_pulse_xy` across all six method/quadrature choices (`expv-pwc`, `expv-pwl`, `expm-pwc`, `expm-pwl`, `evol-pwc`, `evol-pwl`). The regression target states that returned propagators must reproduce two-sided density evolution.

## Behaviour

- Builds a quiet single-proton (`1H`) spin-half system with `test_spin_system` under the `zeeman-hilb` formalism, `zeeman-hilb` approximation `none`, zero Zeeman scalar coupling, and zero magnet.
- Defines noncommuting Pauli control generators (`ops.x`, `ops.y`), a drift Hamiltonian `2*pi*(35*ops.z+17*eye(2))`, and two positive density matrices `rho` and `rho_other`.
- For each of two Hilbert step algorithms (`spin_system.tols.small_matrix` set to `100` and `0`), runs three pulse fixtures: noncommuting shaped pulses, constant pulses, and zero-duration pulses. Piecewise-constant (`pwc`) methods use three amplitude/duration slices; piecewise-linear (`pwl`) methods use four amplitudes with three durations.
- Constructs ordered dense-exponential references via `expm(-1i*H*slice_durs(k))`, with `pwl` references using a Magnus-style midpoint plus a commutator correction `(1i*slice_durs(k)/6)*(H*right-right*H)`.
- For each method/fixture/cutoff combination, checks with tolerances `1e-9` that the returned propagator `P` matches the ordered one-sided product, the returned density matches the two-sided reference `ref_prop*rho*ref_prop'`, `P*rho*P'` reproduces the returned density, `P'*P` is the identity (unitarity), the density trace is 1, the density remains Hermitian, and every trajectory point retains the two-sided state action.
- Verifies output-count invariance: one-output, two-output, and three-output calls return the same density state (exact comparison, tolerance 0), and the propagator is independent of the initial density (`P` equals `Q` computed from `rho_other`, tolerance 0) and reusable on a different initial density (`P*rho_other*P'`).
- Runs one-sided formalism controls for `zeeman-wavef` and `zeeman-liouv` with a nonzero constant pulse over two slices, comparing the returned propagator against `expm(full(-1i*H*sum(slice_durs)))` and checking state-vector reuse `P*initial` with tolerances `1e-8`.
- GPU availability is probed via `exist('canUseGPU','file')` and `canUseGPU()`; when no usable GPU is present, a SKIP message is emitted for the dimension-512 GPU checks.
- Probes sparse `gpuArray` scalar division (`gpuArray(speye(2))/2`); if MATLAB reports that sparse gpuArray matrices are not supported, sparse GPU methods other than `expm-pwc` and `evol-pwc` are skipped with a message, while dense GPU coverage remains complete. Other probe failures are rethrown.
- Lifts the noncommuting spin-half fixture above the GPU dispatch threshold by Kronecker-extending to dimension 512 (`kron(speye(256),...)`), iterating over GPU on/off and sparse/full storage (`spin_system.tols.dense_matrix` set to `1-dense`), with `spin_system.sys.enable` set to `{'gpu'}` or `{}`.
- For the dimension-512 cases, builds independent block-exponential references over two pulse slices and checks the propagator, state, reuse, unitarity against `speye(512)`, trace, Hermiticity, and all trajectory points with tolerances `1e-9`; additionally checks that all returned arrays are gathered to host memory (not `gpuArray`) and that one-output and three-output calls agree with tolerances `1e-12`.

## Inputs and outputs

```matlab
result = test_hilb_pulse_prop()
```

- **Output:** `result` — regression check accumulator created by `new_test_result('kernel/hilb_pulse_prop', ...)`, extended by `test_close` and `test_true` calls, with `result.messages` carrying SKIP notices for unavailable GPU coverage.
- **Input:** none.

## References

- `shaped_pulse_xy` — function under test.
- `new_test_result`, `test_close`, `test_true` — test harness utilities.
- `test_spin_system` — spin system constructor used by the test.
- `pauli`, `operator` — operator generators.
- Source file (complete implementation): https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_hilb_pulse_prop.m
