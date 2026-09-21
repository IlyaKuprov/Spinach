# experiments/nmr_protein/hnco.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/experiments/nmr_protein/hnco.m`
- Signature: `fid=hnco(spin_system,parameters,H,R,K)`
- Total lines: 254

## Purpose

Phase-sensitive HNCO pulse sequence from using the bidirectional propagation method described in The sequence is hard-wired to work on 1H,13C,15N proteins and uses PDB labels to select spins that will be affected by otherwise ideal pulses. F1 is N, F2 is CO, F3 is H. Syntax: fid=hnco(spin_system,parameters,H,R,K)

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
- fid.neg_pos, fid.neg_neg that are used in the subsequ-
- ent States quadrature processing
- Note: spin labels must be set to PDB atom IDs ('CA', 'HA', etc.) in
- sys.labels for this sequence to work properly.

## Implementation structure

- Phase-sensitive HNCO pulse sequence from
- using the bidirectional propagation method described in
- The sequence is hard-wired to work on 1H,13C,15N proteins and uses
- PDB labels to select spins that will be affected by otherwise ideal
- pulses. F1 is N, F2 is CO, F3 is H. Syntax:
- fid=hnco(spin_system,parameters,H,R,K)
- parameters.npoints -a vector of three integers giving the
- number of points in the three temporal
- dimensions, ordered as [t1 t2 t3].
- parameters.sweep -a vector of three real numbers giving
- the sweep widths in the three frequen-
- cy dimensions, ordered as [f1 f2 f3].

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `strcmp()`, `state()`, `operator()`, `step()`, `evolution()`, `coherence()`, `decouple()`, `report()`, `stitch()`, `fieldnames()`, `ismember()`, `ismatrix()`, `all()`, `isfield()`, `isvector()`.
