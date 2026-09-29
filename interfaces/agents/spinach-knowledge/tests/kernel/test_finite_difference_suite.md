# tests/kernel/test_finite_difference_suite.m

**Source:** [https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_finite_difference_suite.m](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_finite_difference_suite.m)

## Purpose

Regression test suite for the finite-difference and spectral differentiation helpers in Spinach. The suite verifies finite-difference weights, finite-difference matrices, Fourier differentiation, Laplacians, FFT differentiation kernels, pseudomodulation, and matrix-exponential directional derivatives against exact simple cases.

## Behaviour

The function announces the test target with `fprintf`, initialises a test result object via `new_test_result` under the identifier `kernel/finite_difference_suite`, and then runs a sequence of checks using `test_close` and `test_true`:

- **`fdweights`**: three-point centred finite-difference weights at zero are checked for orders 0, 1 and 2. The zeroth-order weights must interpolate the centre point (`[0 1 0]`), the first derivative weights must be `[-1/2, 0, 1/2]`, and the second derivative weights must be `[1, -2, 1]`, all with tolerances `1e-14`.
- **`fdmat` and `fdvec`**: a five-point wall finite-difference matrix built with `fdmat(7,5,1,'wall')` must differentiate `x.^2` exactly on a unit grid of seven points, giving `2*x` with tolerances `1e-12`. `fdvec(f,5,1)` must apply the same derivative to a vector.
- **`sgolaydiff`**: local cubic least-squares differentiation must recover the exact derivative of a cubic polynomial `2*x.^3 - x.^2 + 3*x - 1` on the grid `(-5:5)'`, with tolerances `1e-12`. When given a two-column matrix `[f 2*f]`, the columns must be processed independently. The function must reject row-vector signals (error message containing `sample rows`) and even window lengths (error message containing `npoints must be an odd integer`); both rejections are verified with `try/catch` blocks. On a 101-point `2*pi`-periodic grid with a clean sine plus deterministic high-frequency perturbations (`0.03*sin(37*grid) + 0.02*cos(29*grid)`), the Savitzky-Golay derivative error against `cos(grid)` must be smaller than the raw `fdvec` derivative error.
- **`fdmat` (periodic) and `fdlap`**: the periodic first-derivative matrix `fdmat(8,5,1,'pbc')` must annihilate a constant vector (tolerance `1e-14`), and the periodic finite-difference Laplacian `fdlap([5 4],[1 2],3)` must annihilate constants in multiple dimensions (tolerance `1e-12`).
- **`fdkup`**: with an isotropic tensor `eye(3)`, `fdkup([4 4 4],[4 4 4],eye(3),3)` must equal minus one third of the finite-difference Laplacian `fdlap([4 4 4],[4 4 4],3)`, with tolerances `1e-13`.
- **`fourdif` and `fourlap`**: for `N=16`, the first and second Fourier differentiation matrices must exactly differentiate `sin(grid)` to `cos(grid)` and `-sin(grid)` respectively (tolerances `1e-12`). The Fourier Laplacian `fourlap(N,2*pi)` must satisfy `Lf*s = -s`, since `sin(x)` is a Laplacian eigenfunction with eigenvalue `-1` on a `2*pi` periodic interval.
- **`fftdiff`**: the FFT differentiation kernel `fftdiff(1,N,2*pi/N)` must reproduce the same spectral derivative, verified by comparing `real(ifft(fft(s).*kern))` with `cos(grid)` (tolerances `1e-12`).
- **`pseudomodulation`**: the zeroth harmonic with zero amplitude must return the input spectrum unchanged (tolerances `1e-14`). With a small amplitude `pm_amp = 1e-4`, the first harmonic must match the Hyde derivative limit `pm_amp*cos(2*pm_field)` (tolerances `1e-9`/`1e-12`), and the second harmonic must match `(pm_amp^2/4)*sin(2*pm_field)` (tolerances `1e-13`/`1e-14`). The function must reject spectra oriented across columns (error message containing `same number of rows`) and row-vector field axes (error message containing `column vector`); both rejections are verified with `try/catch` blocks.
- **`dirdiff`**: for a spin system with `output='hush'`, empty enable/disable lists, `small_matrix=10` and `prop_chop=1e-14`, and for commuting diagonal matrices `A=diag([1 2])`, `B=diag([3 4])` with `Tstep=0.125`, the zeroth directional derivative must equal the unperturbed propagator `expm(-1i*A*Tstep)` and the first derivative must equal `(-1i*Tstep)*B*P`, both with tolerances `1e-14`.

## Inputs and outputs

```matlab
result=test_finite_difference_suite()
```

**Outputs**

- `result` — regression test result object with explanatory messages, accumulated through the `test_close` and `test_true` checks described above.

**Inputs**

None.

## References

- Source file: [https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_finite_difference_suite.m](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_finite_difference_suite.m)
