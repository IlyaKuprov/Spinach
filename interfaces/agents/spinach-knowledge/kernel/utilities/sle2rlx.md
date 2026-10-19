# kernel/utilities/sle2rlx.m

**Source:** [https://github.com/IlyaKuprov/Spinach/blob/main/kernel/utilities/sle2rlx.m](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/utilities/sle2rlx.m)

## Purpose

Extracts a Redfield-type relaxation superoperator from a stochastic Liouville equation (SLE) generator by Markovian back-projection onto the spin Liouville space. The generator acts on the direct product of a lattice space (Wigner functions of molecular orientation in `gridfree.m`, or chemical species) and the spin Liouville space; the lattice is projected out and its dynamics appears in the spin space as a time-independent relaxation superoperator.

## Behaviour

- Syntax: `[R,L0]=sle2rlx(spin_system,L,pi_vec,one_vec)`.
- `L` is the SLE generator `H+1i*R+1i*K` in the `kron(lattice,spin)` index order of `gridfree.m`, sparse or polyadic (polyadics are inflated). `pi_vec` is the lattice equilibrium distribution and `one_vec` the lattice integration functional, both column vectors; the contraction is `one_vec'*(...)` and the normalisation `one_vec'*pi_vec=1` is applied internally. For `gridfree.m` generators both are the unit vector of the rank-zero Wigner function, which is the first lattice basis element; for chemical exchange they are the equilibrium populations and a vector of ones.
- With the oblique projector `P=|pi><one|(x)1` and `Q=1-P`, Nakajima-Zwanzig partitioning gives the exact spin-space resolvent `[s-L0-C*(s-A)^(-1)*B]^(-1)` with `A=QLQ`, `B=QL|pi>`, `C=<one|LQ`. The Markov limit is the Sylvester equation `A*X-X*L0=-B`, after which `1i*R=C*X`.
- `L0=<one|L|pi>` is the lattice-averaged Liouvillian (the isotropic Hamiltonian for `gridfree.m` generators, the population-weighted average for exchange). Its complex Schur form is taken, and the columns of `X` are obtained in sequence by sparse direct solves of bordered systems `[L-T(n,n)*1, -P_inj; P_con, 0]`, which enforce a zero lattice-equilibrium component of `X` without an explicit basis of the lattice-excitation space, so dense `pi_vec` and `one_vec` are handled at the same cost as sparse ones.
- `R` is returned in the same frame and representation as `L` (no frame transformation is performed; a generator built under rotating-frame assumptions gives a rotating-frame `R`, and comparison with `relaxation.m` needs a lab-frame generator) in the convention of `relaxation.m` (enters the Liouvillian as `1i*R`; dynamic frequency shifts are its imaginary part). When the lattice coupling inside `A` is negligible the result is Redfield theory with spectral densities taken from the lattice generator; the full `A` adds the higher-order corrections. The function does not check the Redfield validity condition (lattice relaxation faster than spin relaxation).
- Any relaxation superoperator already present in `L` (for example when the spin system carries `inter.relaxation`) is averaged into `L0` and is not part of `R`; generators meant for comparison with Redfield theory must be built from a spin system without a relaxation theory.
- Reports lattice and spin dimensions and the 1-norm of `R` through `report`.

## Related

- `kernel/contexts/gridfree.m` builds the SLE generators this function consumes.
- `kernel/relaxation.m` computes the Redfield superoperator that the back-projection reproduces in the fast-motion limit.
- `examples/relaxation_theory/sle2rlx_vs_redfield.m` compares the two for a nitroxide radical and a 15N-1H spin pair.
- `tests/kernel/test_sle2rlx.m` checks two-site exchange against the exact spin line and a CSA-relaxed spin against Redfield theory.
