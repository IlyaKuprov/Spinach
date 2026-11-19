# noise_correlation

Syntax: `[A,aux]=noise_correlation(spin_system,parameters,H,R,K)`

The retained subsystem supplies physical bath IDs. `parameters.noise` explicitly declares Hermitian support/operator terms, `ordering` (`ordered` or `symmetrized`), logical `connected`, and `manifold` (`bath_only`, `alpha`, or `beta`). Times are elapsed seconds; outputs have squared operator units. `aux.moment`, `aux.mean_time`, `aux.mean_zero`, validity, and absolute-error bounds are retained even for connected conditional output. The context centres with weighted ensemble means before additive inclusion-exclusion, including the empty reference. Bath-only GKSL jumps are supported; observed-spin jumps are rejected. Fourier transformation remains separate, with the positive-sign convention in cce_spectrum.

Sparse graph selection closes the support union for unconnected correlations and sampled connected correlations, retaining disconnected operator cross terms. Fixed alpha/beta noise manifolds must be invariant under the retained Hamiltonian.

## Array-free fluctuation descriptors

Instead of noise.terms, supply a nonempty structure array noise.descriptors.
Every element has operators (Lx, Ly, Lz, or E), spins (a nonempty row of
unique physical bath IDs), and coefficient (a finite real scalar in the
fluctuation units). It denotes the coefficient times the sum of those
one-spin operators. The experiment constructs only retained-cluster matrices;
selected exterior spins contribute zero in that retained cluster. All physical
IDs are validated before restriction, including for the empty cluster.
Descriptors and explicit matrix terms are mutually exclusive. Unknown fields,
non-Hermitian labels, observed IDs, duplicate IDs, and non-real coefficients
are rejected. Connected ensemble and Fourier conventions are unchanged.

With no bath dissipator, noise_correlation evaluates the explicitly ordered or symmetrised correlation in Hilbert space using the Hermitian bath eigensystem. It transports the existing unitary round-off estimate. With bath dissipation the source and detector are contracted in Liouville space. The floating-point estimate uses their Frobenius scales and the gain of the computed propagator; source formation and detection still cost arithmetic at zero time. Connected estimates additionally allow for both first moments and their subtraction. These are numerical estimates in squared fluctuation units, not certified error bounds. Optional parameters.spectrum is evaluated after the completed native ensemble and returned in diagnostics.spectrum, not per cluster.

Explicit operator Hermiticity is checked relative to its own Frobenius norm: norm(B-B','fro') must not exceed 100*eps*norm(B,'fro'). There is no absolute unit-scale floor; rescaling a non-Hermitian operator cannot make it acceptable.

Caller-supplied parameters.sequence is rejected, including empty arrays. Use finite_drive, refocusing_train, or the appropriate named echo experiment; internal events are constructed from native physical declarations, never read from that removed field.
