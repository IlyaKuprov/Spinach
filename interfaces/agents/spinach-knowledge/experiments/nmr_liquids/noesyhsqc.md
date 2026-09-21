# experiments/nmr_liquids/noesyhsqc.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/experiments/nmr_liquids/noesyhsqc.m`
- Signature: `fid=noesyhsqc(spin_system,parameters,H,R,K)`
- Total lines: 198

## Purpose

Phase-sensitive NOESY-HSQC sequence described in: The sequence is hard-wired to {F1,F2,F3}={1H,15N,1H} with carbon decoupled throughout. Syntax: fid=noesyhsqc(spin_system,parameters,H,R,K)

## Physical / mathematical content

- Liquid-state pulse sequence implementations. These are production experiment kernels that carry out coherence transfer, mixing, refocusing, decoupling, and indirect evolution on precomputed Hamiltonian/relaxation/kinetics operators.
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
- parameters.J -J-coupling for the HSQC stage magneti-
- sation transfer, Hz.
- parameters.tmix -NOESY stage mixing time, seconds.
- H -Hamiltonian matrix, received from context function
- R -relaxation superoperator, received from context function
- K -kinetics superoperator, received from context function

## Outputs

- fid -a structure with four fields: fid.pos_pos, fid.pos_neg,
- fid.neg_pos, fid.neg_neg that are used in the subsequ-
- ent States quadrature processing
- Notes: the sequence starts with a pure Lz on protons at the mo-
- ment and assumes that the relaxation superoperator is not
- thermalised -the relaxation destination is the zero state.
- The first citation is the NOESY-HMQC precursor to this
- NOESY-HSQC implementation.

## Implementation structure

- Phase-sensitive NOESY-HSQC sequence described in:
- The sequence is hard-wired to {F1,F2,F3}={1H,15N,1H} with carbon
- decoupled throughout. Syntax:
- fid=noesyhsqc(spin_system,parameters,H,R,K)
- parameters.npoints -a vector of three integers giving the
- number of points in the three temporal
- dimensions, ordered as [t1 t2 t3].
- parameters.sweep -a vector of three real numbers giving
- the sweep widths in the three frequen-
- cy dimensions, ordered as [f1 f2 f3].
- parameters.J -J-coupling for the HSQC stage magneti-
- sation transfer, Hz.

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `ismember()`, `decouple()`, `state()`, `operator()`, `step()`, `evolution()`, `coherence()`, `homospoil()`, `report()`, `stitch()`, `fieldnames()`, `ismatrix()`, `all()`, `isfield()`, `isvector()`.
