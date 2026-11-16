# noise_correlation

Syntax: `[A,aux]=noise_correlation(spin_system,parameters,H,R,K)`

The retained subsystem supplies physical bath IDs. `parameters.noise` explicitly declares Hermitian support/operator terms, `ordering` (`ordered` or `symmetrized`), logical `connected`, and `manifold` (`bath_only`, `alpha`, or `beta`). Times are elapsed seconds; outputs have squared operator units. `aux.moment`, `aux.mean_time`, `aux.mean_zero`, validity and absolute-error bounds are retained even for connected conditional output. The context centres with weighted ensemble means before additive inclusion-exclusion, including the empty reference. Bath-only GKSL jumps are supported; observed-spin jumps are rejected. Fourier transformation remains separate, with the positive-sign convention in cce_spectrum.

Sparse graph selection closes the support union for unconnected correlations and sampled connected correlations, retaining disconnected operator cross terms. Fixed alpha/beta noise manifolds must be invariant under the retained Hamiltonian.
