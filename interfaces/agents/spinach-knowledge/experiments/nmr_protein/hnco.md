# experiments/nmr_protein/hnco.m

- Signature: `fid=hnco(spin_system,parameters,H,R,K)`

## Purpose

Phase-sensitive HNCO pulse sequence from [the reported experiment](http://dx.doi.org/10.1016/0022-2364(90)90333-5), using the bidirectional propagation method described in [the 2014 paper](http://dx.doi.org/10.1016/j.jmr.2014.04.002).

## Physical / mathematical content

This sequence is hard-wired for 1H, 13C, and 15N proteins. F1 is N, F2 is CO, and F3 is H. PDB atom labels in `sys.labels` (for example, `CA` and `HA`) select spins affected by otherwise ideal pulses.

## Numerical / algorithmic content

The implementation forms forward and backward propagation stacks, then calls `stitch` to combine them into the four phase/sign combinations returned in `fid`.

## Parameters / inputs

- parameters.npoints -a vector of three integers giving the
- number of points in the three temporal
- dimensions, ordered as [t1 t2 t3].
- parameters.sweep -a vector of three real numbers giving
- the sweep widths in the three frequen-
- cy dimensions, ordered as [f1 f2 f3].
- parameters.tau -the three delays required for the ope-
- ration of the sequence (see the paper)
- in seconds. Reasonable values are
- [2.25e-3, 14e-3, 4e-3]
- parameters.f1_decouple -logical switch controlling proton de-
- coupling during the T1 period.
- H -Hamiltonian matrix, received from context function
- R -relaxation superoperator, received from context function
- K -kinetics superoperator, received from context function

## Outputs

- fid -a structure with four fields: fid.pos_pos, fid.pos_neg,
- fid.neg_pos, fid.neg_neg that are used in the subsequ-- ent States quadrature processing
- Note: spin labels must be set to PDB atom IDs ('CA', 'HA', etc.) in
- sys.labels for this sequence to work properly.
