# kernel/derivatives/fourlap_poly.m

[MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/derivatives/fourlap_poly.m)

`L=fourlap_poly(spin_system,npoints,extents)` is a switch-aware entry point for the periodic Fourier Laplacian. Grid sizes and periods are matching row vectors for one to three coordinates. The operator acts on column-wise vectorisation of an array ordered X, Y, Z.

With polyadics enabled, each coordinate contributes a three-factor FFT second derivative and identities on the other coordinates. Unequal periods and even-grid Nyquist second derivatives are retained. The output is an implicit polyadic that cannot be inflated. Upload it once for GPU actions and use exponential-action propagation for diffusion. Without polyadics, the function returns the existing explicit `fourlap` matrix. The legacy two-argument `fourlap` API remains unchanged.
