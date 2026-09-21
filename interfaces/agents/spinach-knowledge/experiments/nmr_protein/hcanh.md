# experiments/nmr_protein/hcanh.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/experiments/nmr_protein/hcanh.m`
- Signature: `fid=hcanh(spin_system,parameters,H,R,K)`
- Total lines: 272

## Purpose

Protein-specific H(CA)NH experiment (Figure 7.37 of "Protein NMR Spectroscopy", 2nd edition) using pre-set values of J-couplings used in the magnetisation transfer stages. The simulation uses the bidirectional propagation method described in The sequence is hard-wired to work on 1H,13C,15N proteins and uses PDB labels to select spins that will be affected by otherwise ideal pulses. F1 is 1H, F2 is 15N, F3 is 1H. Synt

## Physical / mathematical content

- Protein triple-resonance sequence implementations. They orchestrate heteronuclear coherence transfers across biomolecular spin networks while preserving phase and acquisition conventions.
- Orientation or trajectory averaging is performed numerically, so grid design, weights, and integration error control matter directly to accuracy and runtime.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- Numerical integration over angles or geometry is part of the implementation, so point placement and weights are as important as the local Hamiltonian calculations.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Parameters / inputs

- parameters.npoints -a vector of three integers giving the
- number of points in the three temporal
- dimensions, ordered as [t1 t2 t3].
- parameters.sweep -a vector of three real numbers giving
- the sweep widths in the three frequen-
- cy dimensions, ordered as [f1 f2 f3].
- parameters.spins -isotopes affected by ideal broadband
- pulses, specified as a cell array of
- strings, e.g. {'1H','15N','1H'}.
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

- Protein-specific H(CA)NH experiment (Figure 7.37 of "Protein NMR
- Spectroscopy", 2nd edition) using pre-set values of J-couplings
- used in the magnetisation transfer stages. The simulation uses
- the bidirectional propagation method described in
- The sequence is hard-wired to work on 1H,13C,15N proteins and uses
- PDB labels to select spins that will be affected by otherwise ideal
- pulses. F1 is 1H, F2 is 15N, F3 is 1H. Syntax:
- fid=hcanh(spin_system,parameters,H,R,K)
- parameters.npoints -a vector of three integers giving the
- number of points in the three temporal
- dimensions, ordered as [t1 t2 t3].
- parameters.sweep -a vector of three real numbers giving

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `isfield()`, `ismember()`, `state()`, `operator()`, `strcmp()`, `step()`, `coherence()`, `evolution()`, `decouple()`, `stitch()`, `fieldnames()`, `ismatrix()`, `all()`, `isvector()`, `any()`.
