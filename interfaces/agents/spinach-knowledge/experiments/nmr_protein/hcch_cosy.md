# experiments/nmr_protein/hcch_cosy.m

- Signature: `fid=hcch_cosy(spin_system,parameters,H,R,K)`

## Purpose

HCCH-COSY pulse sequence from Figure 7.26a of the second edition of *Protein NMR Spectroscopy*. It uses the bidirectional propagation method described in [the cited paper](http://dx.doi.org/10.1016/j.jmr.2014.04.002), with ideal pulses selected using PDB atom labels. The sequence is hard-wired for 1H and 13C proteins; F1, F2, and F3 are 1H, 13C, and 1H, respectively.

## Physical / mathematical content

The sequence starts from 1H longitudinal magnetisation, forms positive and negative 1H coherence for F1, and propagates the two branches through the 1H-13C transfer periods. Detection is on 1H; the backward half uses an adjoint propagation and the requested F3 decoupling. The two halves are stitched for the 13C States quadrature and the three dimensions are returned in F3-F2-F1 order.

The code derives `tau_ch = 1/(4 J_ch)` and `tau_cc = 1/(8 J_cc)`, and sets `DELTA = tau_cc - delta`. Evolution uses `L = H + iR + iK`.

## Parameters / inputs

- `parameters.npoints`: three positive integer point counts ordered as [t1 t2 t3].
- `parameters.sweep`: three positive sweep widths in Hz ordered as [f1 f2 f3].
- `parameters.spins`: must be `{'1H','13C','1H'}`.
- `parameters.J_cc`: positive 13C-13C coupling in Hz; the source comment gives 35 Hz as typical.
- `parameters.J_ch`: positive 1H-13C coupling in Hz; the source comment gives 140 Hz as typical.
- `parameters.delta`: positive evolution delay in seconds, shorter than `1/(8 J_cc)`; the source comment gives 1.1 ms as typical.
- `parameters.decouple_f3`: cell array of isotope strings to decouple during F3 acquisition; the source comment gives `{'13C'}` as typical.
- `H`: Hamiltonian matrix; `R`: relaxation superoperator; `K`: kinetics superoperator, supplied by the context function with matching dimensions.
- The spin-system labels must use PDB atom IDs such as CA, HA, and C for selective pulse operations.

## Outputs

Returns `fid.pos_pos`, `fid.pos_neg`, `fid.neg_pos`, and `fid.neg_neg`, the four sign combinations used in subsequent States quadrature processing.

## References

- *Protein NMR Spectroscopy*, 2nd edition, Figure 7.26a.
- [Bidirectional propagation method](http://dx.doi.org/10.1016/j.jmr.2014.04.002)
- [Spin Dynamics Wiki: hcch_cosy.m](https://spindynamics.org/wiki/index.php?title=hcch_cosy.m)
