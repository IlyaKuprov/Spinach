# Periodic DDF regression

The registered suite `kernel/ddf_imaging` tests signed complex signals, not just echo magnitudes. Its single-spin references isolate sample shape, grating orientation, winding order, and diffusion. For the coupled pair it compares the imaging trajectory against a Pauli-matrix Fourier–Galerkin calculation with non-circular mode convolution. Spatial, temporal, and reference-bandwidth refinements distinguish numerical convergence from agreement at one discretisation.
