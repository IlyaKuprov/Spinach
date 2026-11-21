# examples/relaxation_theory/sle2rlx_vs_redfield.m

## Purpose

Compares the relaxation superoperator extracted from a stochastic Liouville equation generator by `sle2rlx` with the Redfield superoperator from `relaxation` for two spin systems: the nitroxide radical of `sle_esr_nitroxide_1.m` (3.5 T, rotational correlation time 5e-11 s) and the 15N-1H pair with dipolar and CSA interactions of `sle_nmr_dd_csa.m` (14.1 T, 5e-9 s).

## Behaviour

- The SLE generator is captured from `gridfree` on a spin system without any relaxation theory, so that the generator contains only the Hamiltonian and the rotational diffusion operator; a second `create` call with the Redfield settings below provides the reference superoperator on the same basis.
- Each reference spin system is built with `inter.relaxation={'redfield'}`, `inter.rlx_keep='labframe'`, `inter.rlx_dfs='keep'`, and `inter.equilibrium='zero'`, so that the Redfield superoperator is the full lab-frame non-secular one with dynamic frequency shifts.
- The SLE generator is obtained from `gridfree` with a capture handle in place of the pulse sequence, with `'labframe'` assumptions so that both calculations see the same lab-frame Hamiltonian, and `parameters.max_rank=4` (rank two is the minimum lattice rank for a rank-two interaction; the nitroxide numbers change by less than 1e-4 relative between ranks two and four).
- The lattice equilibrium vector is the rank-zero Wigner function, the first lattice basis element, as in the `gridfree` state projector.
- Prints normalised matrix elements of both superoperators: longitudinal and transverse rates of the electron and of 15N, the g-hyperfine cross-correlation element between electron coherence and electron-nitrogen product coherence, the 15N-1H cross-relaxation rate, and the CSA-dipolar cross-correlation element between 15N coherence and the 15N-1H antiphase coherence; then the relative Frobenius norm of the difference.
- The residual differences (a few tenths of a percent for the nitroxide, about 1e-4 for the spin pair) are the higher-order corrections to Redfield theory that the full SLE memory kernel contains, not numerical error.
