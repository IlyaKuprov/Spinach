# experiments/nmr_protein/hncaco.m

- Signature: `fid=hncaco(spin_system,parameters,H,R,K)`

## Purpose

Protein-specific HN(CA)CO experiment, Figure 7.41 in the second edition of *Protein NMR Spectroscopy*. It uses preset J-couplings for magnetisation transfer and the bidirectional propagation method described in [the cited paper](http://dx.doi.org/10.1016/j.jmr.2014.04.002). The sequence is hard-wired for 1H, 13C, and 15N proteins; F1, F2, and F3 are 15N, 13C, and 1H, respectively.

## Physical / mathematical content

The forward half starts from positive and negative 15N coherence (the source uses these states to emulate the INEPT block), evolves the indirect 15N dimension with proton decoupling, and applies the CA/CO transfer pulses. The backward half propagates the 1H detection state under the adjoint Liouvillian with 15N decoupling. Stitching the halves yields the four sign combinations for States quadrature; the output dimensions are reordered to F3-F2-F1.

The code derives `tau = 1/(4 J_nh)` and `delta = 1/(2 J_nh)`. `parameters.T` is the indirect 15N evolution delay and `parameters.delta2` is the coherence-transfer delay. Evolution uses `L = H + iR + iK`.

## Parameters / inputs

- `parameters.npoints`: three positive integer point counts ordered as [t1 t2 t3].
- `parameters.sweep`: three positive sweep widths in Hz ordered as [f1 f2 f3].
- `parameters.J_nh`: positive 1H-15N coupling in Hz.
- `parameters.T`: positive evolution delay in seconds; it must be longer than `1/J_nh`.
- `parameters.delta2`: positive coherence-transfer delay in seconds.
- `H`: Hamiltonian matrix; `R`: relaxation superoperator; `K`: kinetics superoperator, supplied by the context function with matching dimensions.
- The spin-system labels must use PDB atom IDs such as CA and C for the sequence's selective operations.

## Outputs

Returns a structure with `fid.pos_pos`, `fid.pos_neg`, `fid.neg_pos`, and `fid.neg_neg`, the four sign combinations used in subsequent States quadrature processing.

## References

- *Protein NMR Spectroscopy*, 2nd edition, Figure 7.41.
- [Bidirectional propagation method](http://dx.doi.org/10.1016/j.jmr.2014.04.002)
- [Spin Dynamics Wiki: hncaco.m](https://spindynamics.org/wiki/index.php?title=hncaco.m)
