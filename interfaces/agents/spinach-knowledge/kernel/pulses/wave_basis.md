# kernel/pulses/wave_basis.m

- MATLAB source: https://github.com/IlyaKuprov/Spinach/blob/main/kernel/pulses/wave_basis.m
- Wiki: https://spindynamics.org/wiki/index.php?title=wave_basis.m
- Signature: `basis_waves=wave_basis(basis_type,n_func,n_points)`

## Purpose

Returns sampled basis functions for pulse-waveform expansion. It constructs the requested functions as rows and returns their orthogonalised sampled vectors as columns of `basis_waves`.

## Basis definitions and discretisation

- `sine_waves`: rows `sin(n*x)` for `n=1:n_func`, with `x` sampled by `linspace(-pi,pi,n_points)`.
- `cosine_waves`: rows `cos((n-1)*x)` on that same interval, so the first row is the constant (zero-frequency) function and subsequent rows begin at frequency 1.
- `legendre`: Legendre polynomials of orders `0:n_func-1` sampled on `linspace(-1,1,n_points)`; each sampled row is first normalised by its 2-norm.

After construction the source applies MATLAB's `orth` to the transpose, making the sampled functions orthogonal as vectors under the discrete representation. The source notes that discretisation means the functions are not precisely orthogonal under the continuous standard scalar product, and that orthogonalisation can flip some functions upside-down. If the requested sampled rows are linearly dependent, the returned column count is smaller than `n_func` and the function errors, directing the caller to reduce `n_func`.

## Inputs and output

- `basis_type` - character string: `sine_waves`, `cosine_waves`, or `legendre`
- `n_func` - positive integer number of requested functions
- `n_points` - positive integer number of discretisation points
- `basis_waves` - matrix with the orthogonalised sampled basis functions in columns

These are dimensionless sampled basis functions, not a pulse-file reader or a pulse amplitude/phase generator. The function exposes no filter parameter.
