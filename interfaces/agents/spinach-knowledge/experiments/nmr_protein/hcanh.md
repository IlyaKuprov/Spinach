# experiments/nmr_protein/hcanh.m

[Canonical MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/experiments/nmr_protein/hcanh.m)

- Signature: `fid=hcanh(spin_system,parameters,H,R,K)`

## Purpose

Protein-specific H(CA)NH, Figure 7.37 of the second edition of *Protein NMR Spectroscopy*. The implementation uses the bidirectional-propagation method described in [the cited paper](http://dx.doi.org/10.1016/j.jmr.2014.04.002). It is hard-wired for 1H, 13C, and 15N; F1, F2, and F3 are 1H, 15N, and 1H.

## Inputs and spin labels

- `spin_system` must use the `sphten-liouv` formalism. Set protein atom labels to PDB atom IDs such as `CA`, `HA`, and `H` so the sequence can select its pulse and state sites.
- `parameters.npoints` is a three-integer vector `[n1 n2 n3]` for `[t1 t2 t3]`; `parameters.sweep` is a three-positive-real vector `[f1 f2 f3]`.
- `parameters.spins` is fixed by the source to `{'1H','15N','1H'}`; it is not a free choice of detected isotopes.
- `H`, `R`, and `K` are same-size matrices supplied by the context function (Hamiltonian, relaxation, and kinetics matrices).
- `parameters.rho0` and `parameters.coil` may be supplied. If absent, the code builds `rho0` from longitudinal magnetisation on labels `HA`, `HA1`, `HA2`, or `HA3`, and builds the proton detection state from label `H`.

## Coherence selection and transfer

The source starts from the labelled H-alpha sites, applies a proton pulse, and retains positive and negative 1H coherence for the first States dimension. Its CA/H/N pulse and evolution blocks implement the H(CA)NH transfer; positive and negative 15N coherence are selected for F2, and 1H single-quantum coherence is used for proton detection in F3. This is the source-supported pathway summary; individual transfer/refocusing operations are encoded in the pulse blocks rather than exposed as a user-selectable pathway.

## Timing and output

The J values are hard-coded: `J_ch=140` and `J_nh=92` (Hz). The corresponding `tau1=abs(1/(4*J_ch))` and `delta1=abs(1/(4*J_ch))` (about 1.79 ms); `tau2=abs(1/(4*J_nh))` and `delta3=abs(1/(4*J_nh))` (about 2.72 ms). The other fixed delays are `delta2 = 12.5 ms` and `delta4 = 23.0 ms`. These delays are in seconds in the implementation; the J values are in Hz.

The returned structure has four FIDs: `fid.pos_pos`, `fid.pos_neg`, `fid.neg_pos`, and `fid.neg_neg`, the two sign choices in each of the first two States dimensions. Each array is permuted to `[n3 n2 n1]`, corresponding to `[t3 t2 t1]` after acquisition/detection and stitching.

## References

- *Protein NMR Spectroscopy*, 2nd edition, Figure 7.37.
- [Bidirectional propagation method](http://dx.doi.org/10.1016/j.jmr.2014.04.002).
- [Spin Dynamics Wiki: hcanh.m](https://spindynamics.org/wiki/index.php?title=hcanh.m).
