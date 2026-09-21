# kernel/coherence.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/coherence.m`
- Signature: `rho=coherence(spin_system,rho,spec)`
- Total lines: 186

## Purpose

Coherence order selection function -keeps only the specified orders of coherence in the state vector. This is useful as an analytical re- placement for complicated phase cycles. Syntax: rho=coherence(spin_system,rho,spec)

## Physical / mathematical content

- This file belongs to the `kernel` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.
- The file uses a Fokker-Planck-style enlarged state space in which spatial or orientational coordinates are promoted to extra dimensions and coupled to spin dynamics through differential operators.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Parameters / inputs

- rho -a state vector or a horizontal stack thereof;
- in zeeman-hilb, a density matrix or a horizon-
- tal stack thereof
- spec -a cell array containing the specification of
- which coherences to keep on which spins. For
- example
- {{'13C',[1 -1]},{'1H',-1}}
- keeps the states that have coherence order
- ((1 OR -1 on 13C) AND (-1 on 1H))
- instead of specific spins, it is possible to
- specify 'electrons', 'nuclei', and 'all'

## Outputs

- rho -the state vector with the undesired orders of
- spin correlations zeroed out
- Note: this function requires sphten-liouv, zeeman-liouv, or zeeman-
- hilb formalism; Fokker-Planck direct products are supported
- in the Liouville space formalisms. In zeeman-hilb, the densi-
- ty matrices are stretched into Liouville space, filtered the-
- re, and folded back.

## Implementation structure

- Coherence order selection function -keeps only the specified orders
- of coherence in the state vector. This is useful as an analytical re-
- placement for complicated phase cycles. Syntax:
- rho=coherence(spin_system,rho,spec)
- rho - a state vector or a horizontal stack thereof;
- in zeeman-hilb, a density matrix or a horizon-
- tal stack thereof
- spec - a cell array containing the specification of
- which coherences to keep on which spins. For
- example
- {{'13C',[1 -1]},{'1H',-1}}
- keeps the states that have coherence order

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `strcmp()`, `lin2lm()`, `false()`, `ischar()`, `cellfun()`, `state_mask()`, `ismember()`, `all()`, `rho()`, `report()`, `vector()`, `iscell()`, `isvector()`, `any()`.
