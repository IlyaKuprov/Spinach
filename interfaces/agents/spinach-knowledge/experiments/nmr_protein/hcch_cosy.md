# experiments/nmr_protein/hcch_cosy.m

- Signature: `fid=hcch_cosy(spin_system,parameters,H,R,K)`

## Purpose

HCCH-COSY pulse sequence from Figure 7.26a of Protein NMR Spectroscopy (2nd edition) using the bidirectional propagation method described in The sequence is hard-wired to work on 1H,13C proteins and uses PDB la- bels to select spins that will be affected by otherwise ideal pulses. F1 is 1H, F2 is C, F3 is H. Syntax: fid=hcch_cosy(spin_system,parameters,H,R,K)

## Physical / mathematical content

- Protein triple-resonance sequence implementations. They orchestrate heteronuclear coherence transfers across biomolecular spin networks while preserving phase and acquisition conventions.
- Orientation or trajectory averaging is performed numerically, so grid design, weights, and integration error control matter directly to accuracy and runtime.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- Numerical integration over angles or geometry is part of the implementation, so point placement and weights are as important as the local Hamiltonian calculations.

## Parameters / inputs

- parameters.npoints -a vector of three integers giving the
- number of points in the three temporal
- dimensions, ordered as [t1 t2 t3].
- parameters.sweep -a vector of three real numbers giving
- the sweep widths in the three frequen-
- cy dimensions, ordered as [f1 f2 f3].
- parameters.J_cc -13C-13C J-coupling to be used for mag-
- netisation transfer, typically 35 Hz
- parameters.J_ch -1H-13C J-coupling to be used for mag-
- netisation transfer, typically 140 Hz
- parameters.delta -evolution delay, see the pulse sequence
- diagram, typically 1.1e-3 seconds.
- parameters.decouple_f3 -list of spins to be decoupled during
- the detection period, typically {'13C'}
- H -Hamiltonian matrix, received from context function
- R -relaxation superoperator, received from context function
- K -kinetics superoperator, received from context function

## Outputs

- fid -a structure with four fields: fid.pos_pos, fid.pos_neg,
- fid.neg_pos, fid.neg_neg that are used in the subsequ-
- ent States quadrature processing
- Note: spin labels must be set to PDB atom IDs ('CA', 'HA', etc.) in
- sys.labels for this sequence to work properly.

## Implementation structure

- HCCH-COSY pulse sequence from Figure 7.26a of Protein NMR Spectroscopy
- (2nd edition) using the bidirectional propagation method described in
- The sequence is hard-wired to work on 1H,13C proteins and uses PDB la-
- bels to select spins that will be affected by otherwise ideal pulses.
- F1 is 1H, F2 is C, F3 is H. Syntax:
- fid=hcch_cosy(spin_system,parameters,H,R,K)
- parameters.npoints -a vector of three integers giving the
- number of points in the three temporal
- dimensions, ordered as [t1 t2 t3].
- parameters.sweep -a vector of three real numbers giving
- the sweep widths in the three frequen-
- cy dimensions, ordered as [f1 f2 f3].
