# kernel/derivatives/fourdif_poly.m

[MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/derivatives/fourdif_poly.m)

`D=fourdif_poly(npoints,order,extent)` constructs a periodic spectral derivative as inverse FFT, a numeric diagonal Fourier multiplier, and FFT polyadic factors. The grid spans the specified period length; derivative units are inverse coordinate units raised to `order`. For even grids, odd derivatives annihilate the Nyquist mode, while even derivatives retain it, matching `fourdif`.

The returned operator acts on column vectors and horizontal state stacks. Transform adjoints include MATLAB FFT normalisation. The numeric multiplier is built once on CPU; upload the operator once with `gpuArray` for GPU actions. Implicit transform factors cannot be inflated to a matrix. Use `fourdif` when an explicit matrix is required.
