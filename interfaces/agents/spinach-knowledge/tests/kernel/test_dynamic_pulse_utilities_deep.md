# tests/kernel/test_dynamic_pulse_utilities_deep.m

**Source**: [https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_dynamic_pulse_utilities_deep.m](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_dynamic_pulse_utilities_deep.m)

## Purpose

Regression test for dynamic pulse utility paths in Spinach. The test verifies that pulse utilities produce finite dynamic outputs and match direct small-matrix references, covering gradient-pulse dynamics, heterodyne filtering, RLC response transforms, Bruker pulse-file writing, finite-RF R-sequence compilation, waveform basis variants, and pulse-shape variants.

## Behaviour

The test announces its target, initialises a regression result via `new_test_result` with the identifier `kernel/dynamic_pulse_utilities_deep`, and then runs the following checks:

- **Gradient pulse**: builds a one-proton Liouville-space spin system with a carrier frequency (`local_liouv_system(1.0)`), zero `Lz` operator, and an `Lx` state. `grad_pulse` is compared against a direct small-matrix exponential reference (`grad_pulse_ref`) with tolerance `1e-10`, and a non-zero response check (`norm(rho_grad-rho,2)>1e-8`) confirms that a finite gradient changes transverse magnetisation in the one-spin carrier frame.
- **Gradient sandwich**: `grad_sandw` is compared against `grad_sandw_ref` using a propagator built from a `2*pi*50` Ly operator with a `2.0e-4` duration, gradient amplitudes `[3.0 -2.0]`, slice length `0.10`, durations `[1.0e-4 1.4e-4]`, and factors `[0.7 0.9]`, with tolerance `1e-10` and a non-zero response check against the propagated state.
- **Heterodyne**: mixes a unit cosine at 1000 Hz (4096 points, `dt=1e-4`) with its carrier. Checks that the steady-state in-phase mean is 1 and the quadrature mean is 0 (tolerance `5e-2`, steady window from index 200 to end minus 200), and that all output samples are finite.
- **RLC response transforms**: exercises `restrans` with `omega=2*pi*800`, `Q=8`, `dt=1e-3`, and order 2 for the `'pwc'`, `'pwl'`, and `'pwl_tsc'` models. Verifies finite column waveforms and positive time steps for `'pwc'` and `'pwl'`, and that `'pwl_tsc'` returns the same waveform samples (tolerance `1e-12`) and time step (tolerance `1e-15`) as `'pwl'`.
- **Bruker pulse writing**: writes a temporary Bruker pulse file via `bruker_write` with amplitudes `[1;0;-1]`, phases `[0;1;0]`, and duration `2.5e-6`, then checks that the file exists, declares `##NPOINTS=3`, declares `##$SHAPE_LENGTH=7.5` (total pulse length in microseconds), and is terminated by `##END`. A `onCleanup` handler deletes the temporary file.
- **Finite-RF R-sequence compilation**: uses a Hilbert-space system (`local_hilb_system`) with Pauli operators. For a `'180_pulse'` with phases `[0;pi/2;0]`, amplitude `7.5`, and duration `0.03`, checks the index map `T=[1;2;1]` (tolerance `1e-15`) and that the compiled propagators match `expm(-1i*pulse_amp*Sx*pulse_dur)` and `expm(-1i*pulse_amp*Sy*pulse_dur)` (tolerance `1e-14`). For a `'90270_pulse'` with phases `[0;pi/2]` and durations `[0.01 0.02]`, checks the index map `T=[1;2]` and the corresponding exponentials for each duration.
- **Waveform basis variants**: for `{'sine_waves','cosine_waves','legendre'}`, builds `wave_basis(type,3,17)` and checks that the Gram matrix equals `eye(3)` (tolerance `1e-12`) and that the size is `[17 3]`.
- **Pulse-shape variants**: checks `pulse_shape('rectangular',5)` equals `ones(1,5)`; `pulse_shape('sinc3',5)` equals `pi*sinc(linspace(-3,3,5))`; `pulse_shape('sinc5',5)` equals `pi*sinc(linspace(-5,5,5))`; and `pulse_shape('gaussian',5)` equals `exp(-(time_grid.^2)/2)/(sqrt(2*pi)*sqrt(2))` over `linspace(-2,2,5)`, all with tolerance `1e-15`.

Independent small-matrix auxiliary-exponential references provide the expected single-gradient and gradient-sandwich propagations in both Hilbert and Liouville formalisms.

## Inputs and outputs

**Syntax**

```matlab
result = test_dynamic_pulse_utilities_deep()
```

**Outputs**

- `result` — regression test result with explanatory messages.

The function takes no inputs.

## References

- Edwards, D. H. et al. (1992) `grad_pulse` and `grad_sandw` gradient propagation methods, as referenced in the test messages.
- Spinach pulse utility functions: `grad_pulse`, `grad_sandw`, `heterodyne`, `restrans`, `bruker_write`, `rseq_compiler`, `wave_basis`, `pulse_shape`.
