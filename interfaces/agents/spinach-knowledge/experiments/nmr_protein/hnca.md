# experiments/nmr_protein/hnca.m

- Signature: `fid=hnca(spin_system,parameters,H,R,K)`

## Purpose

Protein-specific HNCA experiment, Figure 7.31a in the second edition of *Protein NMR Spectroscopy*. It uses preset J-couplings for magnetisation transfer and the bidirectional propagation method described in [the cited paper](http://dx.doi.org/10.1016/j.jmr.2014.04.002). The sequence is hard-wired for 1H, 13C, and 15N proteins; F1, F2, and F3 are 15N, 13C, and 1H, respectively.

## Physical / mathematical content

The forward half begins with NH-proton magnetisation, transfers it to 15N, and selects positive and negative 15N coherence for F1 States quadrature. The backward half propagates the 1H detection state under the adjoint Liouvillian, with 15N decoupling during F3. Stitching the halves forms the four sign combinations; the output dimensions are reordered to F3-F2-F1.

The source hard-codes J_NH = 92 Hz and J_NCA = 11.5 Hz, with transfer delays derived from those couplings. Evolution uses `L = H + iR + iK`.

## Parameters / inputs

- `parameters.npoints`: three positive integer point counts ordered as [t1 t2 t3].
- `parameters.sweep`: three positive sweep widths in Hz ordered as [f1 f2 f3].
- `parameters.rho0`: optional initial state; if omitted, the source constructs it from protons labelled H.
- `parameters.coil`: optional detection state; if omitted, the source constructs it from protons labelled H.
- `H`: Hamiltonian matrix; `R`: relaxation superoperator; `K`: kinetics superoperator, supplied by the context function with matching dimensions.
- The spin-system labels must use PDB atom IDs such as CA, HA, and C for the selective pulse operations.

## Outputs

Returns a structure with `fid.pos_pos`, `fid.pos_neg`, `fid.neg_pos`, and `fid.neg_neg`, the four sign combinations used in subsequent States quadrature processing.

## References

- *Protein NMR Spectroscopy*, 2nd edition, Figure 7.31a.
- [Bidirectional propagation method](http://dx.doi.org/10.1016/j.jmr.2014.04.002)
- [Spin Dynamics Wiki: hnca.m](https://spindynamics.org/wiki/index.php?title=hnca.m)
