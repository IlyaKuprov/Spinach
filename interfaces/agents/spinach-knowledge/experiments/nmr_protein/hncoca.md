# experiments/nmr_protein/hncoca.m

- Signature: `fid=hncoca(spin_system,parameters,H,R,K)`

## Purpose

Phase-sensitive HN(CO)CA pulse sequence from [the reported experiment](http://dx.doi.org/10.1007/BF01874573), using the bidirectional propagation method described in [the 2014 paper](http://dx.doi.org/10.1016/j.jmr.2014.04.002).

## Physical / mathematical content

This sequence is hard-wired for 1H, 13C, and 15N proteins. F1 is N, F2 is CA, and F3 is H. PDB atom labels in `sys.labels` select spins affected by otherwise ideal pulses.

## Numerical / algorithmic content

The implementation uses bidirectional propagation and `stitch` calls to combine forward and backward evolution into the phase/sign components of `fid`.

## Parameters / inputs

- parameters.npoints -a vector of three integers giving the
- number of points in the three temporal
- dimensions, ordered as [t1 t2 t3].
- parameters.sweep -a vector of three real numbers giving
- the sweep widths in the three frequen-
- cy dimensions, ordered as [f1 f2 f3].
- parameters.tau -the four delays required for the ope-
- ration of the sequence (see the paper)
- in seconds. Reasonable values are
- [2.25e-3, 2.75e-3, 8.00e-3, 7.00e-3]
- H -Hamiltonian matrix, received from context function
- R -relaxation superoperator, received from context function
- K -kinetics superoperator, received from context function

## Outputs

- fid -a structure with four fields: fid.pos_pos, fid.pos_neg,
- fid.neg_pos, fid.neg_neg that are used in the subsequ-
- ent States quadrature processing
- Note: spin labels must be set to PDB atom IDs ('CA', 'HA', etc.) in- sys.labels for this sequence to work properly.
