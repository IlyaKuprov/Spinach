# experiments/nmr_protein/hcch_cosy.m

[Canonical MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/experiments/nmr_protein/hcch_cosy.m)

- Signature: `fid=hcch_cosy(spin_system,parameters,H,R,K)`

## Purpose

HCCH-COSY, Figure 7.26a of the second edition of *Protein NMR Spectroscopy*. The implementation uses the bidirectional-propagation method described in [the cited paper](http://dx.doi.org/10.1016/j.jmr.2014.04.002) and is hard-wired for 1H and 13C. F1, F2, and F3 are 1H, 13C, and 1H.

## Inputs and spin labels

- `spin_system` must use the `sphten-liouv` formalism. PDB atom labels such as `CA`, `HA`, and `C` identify pulse sites, including the carbonyl-carbon site labelled `C`.
- `parameters.spins={'1H','13C','1H'}` is mandatory; the grumbler rejects a missing or different channel selector.
- `parameters.npoints` is a three-integer vector `[n1 n2 n3]` for `[t1 t2 t3]`; `parameters.sweep` is a three-positive-real vector `[f1 f2 f3]`.
- `parameters.J_cc` and `parameters.J_ch` are the 13C-13C and 1H-13C couplings in Hz. Header examples are 35 Hz and 140 Hz, respectively.
- `parameters.delta` is a positive pulse-sequence evolution delay in seconds (header example `1.1e-3`). It must be strictly less than `1/(8*parameters.J_cc)`; equality or a longer delay is rejected before `DELTA=tau_cc-parameters.delta` can become non-positive.
- `parameters.decouple_f3` lists nuclei to decouple during detection; the header example is `{'13C'}`.
- `H`, `R`, and `K` are same-size matrices supplied by the context function (Hamiltonian, relaxation, and kinetics matrices).

## Coherence selection and transfer

The source begins with 1H longitudinal magnetisation, creates positive and negative 1H coherence for F1, and evolves through the 1H-13C and 13C-13C coupling delays. The pulse blocks include a carbonyl-selective pulse on atoms labelled `C`. For the reverse/detection half, the code selects positive and negative 13C coherence for F2 and uses 1H single-quantum detection for F3. This is the HCCH-COSY correlation encoded by the source; the four States sign combinations are returned separately.

## Timing and output

The source sets `tau_ch=abs(1/(4*parameters.J_ch))`, `tau_cc=abs(1/(8*parameters.J_cc))`, and `DELTA=tau_cc-parameters.delta`. With the header's typical values, these are about 1.79 ms, 3.57 ms, and 2.47 ms, respectively. Couplings are in Hz and delays in seconds. `delta` is passed into the sequence's stitched pulse block and its adjoint reverse evolution.

The returned structure has `fid.pos_pos`, `fid.pos_neg`, `fid.neg_pos`, and `fid.neg_neg`. Each FID is permuted to `[n3 n2 n1]`, or `[t3 t2 t1]` in acquisition order.

## References

- *Protein NMR Spectroscopy*, 2nd edition, Figure 7.26a.
- [Bidirectional propagation method](http://dx.doi.org/10.1016/j.jmr.2014.04.002).
- [Spin Dynamics Wiki: hcch_cosy.m](https://spindynamics.org/wiki/index.php?title=hcch_cosy.m).
