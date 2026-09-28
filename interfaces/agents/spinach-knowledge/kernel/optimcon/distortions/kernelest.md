# kernel/optimcon/distortions/kernelest.m

- Signature: `h=kernelest(x,y,ker_len,method,align,lambda)`

## Purpose

Estimate an FIR convolution kernel from input and output signal samples on the same uniform grid.

## Parameters / inputs

- `x` — numeric vector of input samples.
- `y` — numeric vector of output samples, with the same length as `x`.
- `ker_len` — positive integer kernel length, in taps.
- `method` — solution method: `'backslash'` (default), `'pinv'`, `'svd'`, or `'tikh'`.
- `align` — output alignment: `'causal'` (default) or `'same'`.
- `lambda` — positive real scalar Tikhonov parameter for `'tikh'`; defaults to `1e-6`.

## Output

- `h` — estimated convolution kernel.

## Numerical method

The function builds a Toeplitz convolution matrix from `x`, then selects either its first `numel(x)` rows for `'causal'` alignment or a central block starting at `floor(ker_len/2)+1` for `'same'` alignment. It solves the resulting system with MATLAB backslash, a pseudoinverse, a truncated-SVD pseudoinverse, or Tikhonov regularisation, according to `method`. The truncated-SVD method discards singular values at or below `max(size(sys_mat))*eps(max(s_vals))`; the Tikhonov method solves `(sys_mat'*sys_mat+lambda*eye(ker_len))*h=sys_mat'*y`.

The function checks input types, sample counts, kernel length, and `lambda`, and rejects unknown methods or alignment choices.

## Reference

- [Spinach documentation for `kernelest.m`](https://spindynamics.org/wiki/index.php?title=kernelest.m)