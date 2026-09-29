# kernel/optimcon/distortions/kernelest.m

- Signature: h=kernelest(x,y,ker_len,method,align,lambda)
- MATLAB source: [kernelest.m](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/optimcon/distortions/kernelest.m)

## Purpose

Estimates an FIR convolution kernel from input samples x and output samples y on the same uniform grid. It builds a Toeplitz convolution matrix from x and solves the selected linear system; it does not return a gradient or adjoint.

## Inputs and defaults

- x and y are numeric vectors with equal numbers of samples. The function reshapes each to a column; the header describes samples on a uniform grid. No time step or physical units are supplied by this interface.
- ker_len is a positive integer tap count.
- method may be 'backslash' (default), 'pinv', 'svd', or 'tikh'. The choice is lowercased before dispatch, so letter case does not change a recognised choice.
- align may be 'causal' (default) or 'same', also dispatched case-insensitively.
- lambda defaults to 1e-6 and must be a positive real scalar. The value is used by 'tikh'; the source validates it for every method.

The method and alignment inputs must be character vectors or strings. The source does not explicitly require x or y to be real or finite. Unrecognised method or alignment values raise an error during dispatch.

## Matrix construction and methods

The convolution matrix has ker_len columns and is built from the samples of x, followed by zeros. For 'causal' alignment, the first numel(x) rows are selected. For 'same', a central block of numel(x) rows begins at floor(ker_len/2)+1.

- 'backslash' uses MATLAB's backslash solve (labelled least-squares in the source).
- 'pinv' applies the pseudoinverse, giving a minimum-norm least-squares estimate.
- 'svd' forms an economy SVD and uses the truncated pseudoinverse: singular values are retained only when they exceed max(size(sys_mat))*eps(max(s_vals)).
- 'tikh' solves (sys_mat'*sys_mat+lambda*eye(ker_len))*h=sys_mat'*y, equivalent to a Tikhonov-regularised least-squares fit with squared residual plus lambda times squared kernel norm.

## Reference

- Spinach documentation: https://spindynamics.org/wiki/index.php?title=kernelest.m
