# tests/kernel/test_sle2rlx.m

## Purpose

Regression test for `sle2rlx`, the Markovian back-projection of stochastic Liouville equation generators onto the spin Liouville space.

## Behaviour

- Two-site chemical exchange of a single coherence with unequal populations: the generator is `diag([omega_a omega_b])+1i*K` with the exchange matrix `K`, the lattice vectors are the equilibrium populations and a vector of ones. The lattice-averaged Liouvillian must be the population-weighted frequency; `L0+1i*R` must match the exact eigenvalue of the generator nearest to it to within `(delta_omega/k_ex)^2`; and `real(R)` must match the Swift-Connick rate `-p_a*p_b*delta_omega^2/k_ex` to the same tolerance.
- Single 13C spin with an anisotropic chemical shift at 14.1 T and a correlation time of 1e-11 s: the SLE generator is captured from `gridfree` in the lab frame with `max_rank=4`; the averaged Liouvillian must equal the lab-frame isotropic Hamiltonian, and the back-projected superoperator must match the lab-frame Redfield superoperator with dynamic frequency shifts kept to within 1e-4 of its 1-norm, the size of the second-order corrections.
- Registered in `tests/lib/test_manifest.m` as `kernel/sle2rlx`.
