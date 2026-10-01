# kernel/derivatives/fourdif_fft.m

[MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/derivatives/fourdif_fft.m)

`D=fourdif_fft(npoints)` represents the first angular derivative on a uniform periodic grid spanning `2*pi`. It uses the existing `fftdiff` spectral kernel and applies it through FFTs along the first dimension of full `npoints`-by-`n` blocks. Complex states are retained, not projected to their real parts. On even grids the Nyquist derivative is zero, matching `fourdif(npoints,1)` for complex as well as real states.

The result is a `matfree` forward/adjoint core, except for the one-point grid where it is scalar zero. The first derivative is real and skew-adjoint. Use the core inside a polyadic with an identity on the spin subspace to apply rotor turning without opening its Kronecker matrix. This is the same spectral differentiation discretisation, not a reduction in rotor bandwidth. Actions support CPU and GPU blocks; the core does not construct an explicit derivative matrix.
