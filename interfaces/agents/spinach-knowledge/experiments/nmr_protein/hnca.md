# experiments/nmr_protein/hnca.m

[Canonical MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/experiments/nmr_protein/hnca.m)

- Signature: `fid=hnca(spin_system,parameters,H,R,K)`

## Purpose

Protein-specific HNCA, Figure 7.31a of the second edition of *Protein NMR Spectroscopy*. The implementation uses the bidirectional-propagation method described in [the cited paper](http://dx.doi.org/10.1016/j.jmr.2014.04.002). It is hard-wired for 1H, 13C, and 15N; F1, F2, and F3 are 15N, 13C, and 1H.

## Inputs and spin labels

- `spin_system` must use the `sphten-liouv` formalism. Set PDB atom labels such as `H`, `CA`, and `C` so the sequence can identify the proton, alpha-carbon, and carbonyl sites.
- `parameters.npoints` is a three-integer vector `[n1 n2 n3]` for `[t1 t2 t3]`; `parameters.sweep` is a three-positive-real vector `[f1 f2 f3]`.
- `H`, `R`, and `K` are same-size matrices supplied by the context function (Hamiltonian, relaxation, and kinetics matrices). The source header exposes no `parameters.spins` isotope-cell-array option; the isotope channels are fixed by the sequence.
- `parameters.rho0` and `parameters.coil` may be supplied. If absent, the code builds the initial longitudinal state and proton detection state from sites labelled `H`.

## Coherence selection and transfer

The default initial state is on the labelled amide protons. A 1H pulse and the first transfer block create proton coherence; the sequence then selects positive and negative 15N coherence for F1. It refocuses the first evolution interval with H, CA, and CO pulses, applies the CA-transfer delay, and uses the CA/H pulse block before the reverse half. The code selects positive and negative coherence on the labelled CA carbons for F2 and uses 1H single-quantum detection for F3. This is the HNCA transfer implemented by the source, with CA and CO sites selected by their PDB labels.

## Timing and output

The couplings are fixed in the source: `J_nh=92` and `J_nca=11.5` (Hz). The corresponding delays are `tau=abs(1/(4*J_nh))` (about 2.72 ms) and `delta=abs(1/(4*J_nca))` (about 21.74 ms). Couplings are in Hz and delays in seconds.

The returned structure has `fid.pos_pos`, `fid.pos_neg`, `fid.neg_pos`, and `fid.neg_neg`, the four States sign combinations across F1 and F2. Each FID is permuted to `[n3 n2 n1]`, or `[t3 t2 t1]` in acquisition order.

## References

- *Protein NMR Spectroscopy*, 2nd edition, Figure 7.31a.
- [Bidirectional propagation method](http://dx.doi.org/10.1016/j.jmr.2014.04.002).
- [Spin Dynamics Wiki: hnca.m](https://spindynamics.org/wiki/index.php?title=hnca.m).
