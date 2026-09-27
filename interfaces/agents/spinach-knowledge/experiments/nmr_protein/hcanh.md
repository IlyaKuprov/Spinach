# experiments/nmr_protein/hcanh.m

- Signature: `fid=hcanh(spin_system,parameters,H,R,K)`

## Purpose

Protein-specific H(CA)NH experiment, Figure 7.37 in the second edition of *Protein NMR Spectroscopy*. It uses preset J-couplings for magnetisation transfer and the bidirectional propagation method described in [the cited paper](http://dx.doi.org/10.1016/j.jmr.2014.04.002). The sequence is hard-wired for 1H, 13C, and 15N proteins; F1, F2, and F3 are 1H, 15N, and 1H, respectively.

## Physical / mathematical content

The forward half starts from CA-proton magnetisation, selects positive and negative 1H coherence for F1 States quadrature, and carries out the transfer and refocusing steps. The backward half propagates the 1H detection state under the adjoint Liouvillian; stitching the two halves produces the four sign combinations. The pulse sequence uses ideal broadband pulses selected by PDB atom labels and decouples 13CO during the indicated transfer periods.

The source hard-codes J_CH = 140 Hz and J_NH = 92 Hz, with delays derived from those couplings and additional fixed delays of 12.5 ms and 23.0 ms. Evolution uses `L = H + iR + iK`.

## Parameters / inputs

- `parameters.npoints`: three positive integer point counts ordered as [t1 t2 t3].
- `parameters.sweep`: three positive sweep widths in Hz ordered as [f1 f2 f3].
- `parameters.spins`: required to be `{'1H','15N','1H'}`.
- `parameters.rho0`: optional initial state; if omitted, the source constructs it from protons with labels HA, HA1, HA2, or HA3.
- `parameters.coil`: optional detection state; if omitted, the source constructs it from protons labelled H.
- `H`: Hamiltonian matrix; `R`: relaxation superoperator; `K`: kinetics superoperator, supplied by the context function with matching dimensions.
- The spin-system labels must use PDB atom IDs such as CA, HA, and H for the sequence's selective operations.

## Outputs

Returns a structure with `fid.pos_pos`, `fid.pos_neg`, `fid.neg_pos`, and `fid.neg_neg`, the four sign combinations used in subsequent States quadrature processing.

## References

- *Protein NMR Spectroscopy*, 2nd edition, Figure 7.37.
- [Bidirectional propagation method](http://dx.doi.org/10.1016/j.jmr.2014.04.002)
- [Spin Dynamics Wiki: hcanh.m](https://spindynamics.org/wiki/index.php?title=hcanh.m)
