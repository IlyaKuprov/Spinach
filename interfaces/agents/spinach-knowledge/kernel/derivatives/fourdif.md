# kernel/derivatives/fourdif.m

[Direct MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/derivatives/fourdif.m) · [Spinach Wiki documentation](https://spindynamics.org/wiki/index.php?title=fourdif.m)

## Purpose and signature

`[x,DM] = fourdif(spin_system,N,m)` returns Fourier spectral differentiation data on the canonical periodic grid with `N` equally spaced points in `[0,2*pi)`.

## Inputs and outputs

- `spin_system`: Spinach system whose `sys.enable` option list selects the representation.
- `N`: positive real integer scalar, the number of grid points.
- `m`: positive real integer scalar, the derivative order.
- `x`: `N`-by-1 column of points `x_j = 2*pi*j/N`, for `j=0,...,N-1`.
- `DM`: `N`-by-`N` explicit differentiation matrix or an implicit polyadic, formed only when the caller requests a second output.

For example, `[x,DM]=fourdif(spin_system,8,1)` requests eight grid points and an 8-by-8 first-derivative matrix. The coordinate spans one period; the matrix is defined with respect to this canonical coordinate, not automatically rescaled to a different interval.

## Construction

For first and second derivatives, the source uses explicit first-column formulae, with the parity-dependent cotangent/cosecant forms for the first derivative and the flipping trick cited below for improved accuracy. For orders above two, it uses the discrete-Fourier construction. It obtains the first row and column and assembles the matrix with `toeplitz`.

The source also contains a zeroth-order identity branch, while its input guard accepts only positive `m`; the documented accepted call therefore has `m >= 1`.

## Validation and reference

The guard requires `N` and `m` each to be numeric, finite, real, scalar, and a positive integer. The grid spacing used in the construction is `2*pi/N`.

- S.C. Reddy and J.A.C. Weideman, [doi:10.1137/0916073](http://dx.doi.org/10.1137/0916073).
- [Spinach Wiki documentation](https://spindynamics.org/wiki/index.php?title=fourdif.m)

## Polyadic option

With `polyadic` in `spin_system.sys.enable`, the second output consists of inverse FFT, a numeric CPU diagonal multiplier, and FFT polyadic factors. Transform adjoints include MATLAB FFT normalisation. Even-grid odd derivatives annihilate the Nyquist mode, while even derivatives retain it. The canonical period is still `2*pi`; scale the derivative by `(2*pi/extent)^m` for another positive period. The output acts on column vectors and horizontal state stacks. Upload it once with `gpuArray` for GPU actions; implicit transform factors cannot be inflated. Request only the first output for grid points without constructing either derivative representation.
