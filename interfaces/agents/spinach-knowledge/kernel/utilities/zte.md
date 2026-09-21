# kernel/utilities/zte.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/utilities/zte.m`
- Signature: `projector=zte(spin_system,L,rho,nstates)`
- Total lines: 174

## Purpose

Zero track elimination function. Inspects the first few steps in the system trajectory and drops the states that did not get populated to a user-specified tolerance. Syntax: projector=zte(spin_system,L,rho,nstates)

## Physical / mathematical content

- General mathematical and infrastructure utilities. This area contains finite differences, perturbation theory, graph algorithms, spectral densities, tensor algebra, hash/report helpers, and other reusable numerical components.
- Propagation is accelerated with a Krylov-subspace method, replacing direct matrix exponentiation by projection into a much smaller Arnoldi/Lanczos-type subspace.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- A Krylov-subspace or Arnoldi construction is used to avoid forming or exponentiating very large dense propagators directly.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `size()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Parameters / inputs

- L -the Liouvillian to be used for time
- propagation
- rho -the initial state to be used for
- time propagation
- nstates -if this parameter is specified, only
- nstates most populated states are kept,
- irrespective of the tolerance parameter
- Output:
- projector -projector matrix into the reduced space,
- to be used as follows:
- L_reduced=P'*L*P
- rho_reduced=P'*rho;
- Note: default tolerance may be altered by setting sys.tols.zte_tol
- variable before calling create.m
- Note: further information on how this function works is available
- in IK's JMR paper on the subject
- Note: if tiny interactions or nearly equivalent spins are present,
- it is best to disable zero track elimination by adding 'zte'
- to the sys.disable cell array.

## Implementation structure

- Zero track elimination function. Inspects the first few steps in the
- system trajectory and drops the states that did not get populated to
- a user-specified tolerance. Syntax:
- projector=zte(spin_system,L,rho,nstates)
- L -the Liouvillian to be used for time
- propagation
- rho -the initial state to be used for
- time propagation
- nstates -if this parameter is specified, only
- nstates most populated states are kept,
- irrespective of the tolerance parameter
- Output:

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `ismember()`, `report()`, `nnz()`, `cheap_norm()`, `isinf()`, `exist()`, `num2str()`, `trajectory()`, `step()`, `true()`, `zero_track_mask()`, `index()`, `false()`, `speye()`, `projector()`.
