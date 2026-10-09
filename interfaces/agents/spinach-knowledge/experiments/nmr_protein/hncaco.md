# experiments/nmr_protein/hncaco.m

[Canonical MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/experiments/nmr_protein/hncaco.m)

- Signature: `fid=hncaco(spin_system,parameters,H,R,K)`

## Purpose

Protein-specific HN(CA)CO, Figure 7.41 of the second edition of *Protein NMR Spectroscopy*. The implementation uses the bidirectional-propagation method described in [the cited paper](http://dx.doi.org/10.1016/j.jmr.2014.04.002). It is hard-wired for 1H, 13C, and 15N; F1, F2, and F3 are 15N, 13C, and 1H.

## Inputs and spin labels

- `spin_system` must use the `sphten-liouv` formalism. PDB atom labels such as `CA`, `C`, and `H` identify the alpha-carbon, carbonyl-carbon, and proton sites.
- `parameters.npoints` is a three-integer vector `[n1 n2 n3]` for `[t1 t2 t3]`; `parameters.sweep` is a three-positive-real vector `[f1 f2 f3]`.
- `parameters.J_nh` is the 1H-15N coupling in Hz; `parameters.T` is the indirect 15N evolution delay in seconds; `parameters.delta2` is a coherence-transfer delay in seconds. The code requires positive values and `T > 1/J_nh`.
- `H`, `R`, and `K` are same-size matrices supplied by the context function (Hamiltonian, relaxation, and kinetics matrices). The source header exposes no isotope-cell-array option; the channels are fixed in the sequence.

## Coherence selection and transfer

The code initialises positive and negative 15N states to emulate the INEPT block, uses 15N coherence for F1, and applies CA- and CO-selective pulse/evolution blocks. Its reverse half uses proton detection; positive and negative coherence are selected on the carbonyl carbons labelled `C` for F2. Thus the encoded HN(CA)CO experiment reports the 15N/CO dimensions with 1H detection. The source specifies the pulse operations and site selections but does not annotate a separate coherence order for every internal transfer delay.

## Timing and output

From the supplied `J_nh`, the source calculates `tau=abs(1/(4*parameters.J_nh))` and `delta=abs(1/(2*parameters.J_nh))`; J is in Hz and the resulting times are in seconds. `T` is divided into two `T/2` intervals around the N/CA inversion block. The supplied `delta2` is used for the two coherence-transfer intervals around the CO/CA pulses.

The returned structure has `fid.pos_pos`, `fid.pos_neg`, `fid.neg_pos`, and `fid.neg_neg`, the four States sign combinations across F1 and F2. Each FID is permuted to `[n3 n2 n1]`, or `[t3 t2 t1]` in acquisition order.

## References

- *Protein NMR Spectroscopy*, 2nd edition, Figure 7.41.
- [Bidirectional propagation method](http://dx.doi.org/10.1016/j.jmr.2014.04.002).
- [Spin Dynamics Wiki: hncaco.m](https://spindynamics.org/wiki/index.php?title=hncaco.m).
