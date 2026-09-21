# kernel/states/deut_pair.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/states/deut_pair.m`
- Signature: `[S,T,Q,Tc,Qc]=deut_pair(spin_system,spin_a,spin_b,options)`
- Total lines: 198

## Purpose

All possible states of a spin-1 pair, classified by the total spin into singlet, triplet, and quartet. Syntax: [S,T,Q,Tc,Qc]=deut_pair(spin_system,spin_a,spin_b,options)

## Physical / mathematical content

- State-construction utilities. These routines build equilibrium states, singlets, triplets, partner-state expansions, and physically meaningful density operators in the active basis.
- The relevant state manifold is the singlet/triplet decomposition, where permutation symmetry controls selection rules, relaxation susceptibility, and convertibility to ordinary magnetisation.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Parameters / inputs

- spin_a -the number of the first spin
- spin_b -the number of the second spin
- options.dephasing -set to 1 to eliminate the states that are
- not stationary under the Zeeman Hamilto-
- nian of two inequivalent spins, i.e. to
- keep only zero-projection products; the
- default is to keep everything

## Outputs

- S -singet state density matrix (Hilbert space)
- or state vector (Liouville space)
- T -triplet state density matrices (Hilbert space)
- or state vectors (Liouville space), ordered in
- a cell array as {T+,T0,T-}
- Q -quintet state density matrices (Hilbert space)
- or state vectors (Liouville space), ordered in
- a cell array as {Q++,Q+,Q0,Q-,Q--}
- Tc -coherences between triplet states:
- {T0 -> T-, T+ -> T0, T--> T0, T0 -> T+}
- Qc -coherences between quintet states:
- {Q--> Q--, Q0 -> Q-, Q+ -> Q0, Q++ -> Q+, ...
- Q---> Q-, Q--> Q0, Q0 -> Q+, Q+ -> Q++ }
- WARNING: the states above are NOT irreducible spherical tensors -
- Bargon just kroneckered up some Zeeman states and gave
- them what looked to him like reasonable labels.

## Implementation structure

- All possible states of a spin-1 pair, classified by the total spin
- into singlet, triplet, and quartet. Syntax:
- [S,T,Q,Tc,Qc]=deut_pair(spin_system,spin_a,spin_b,options)
- spin_a -the number of the first spin
- spin_b -the number of the second spin
- options.dephasing -set to 1 to eliminate the states that are
- not stationary under the Zeeman Hamilto-
- nian of two inequivalent spins, i.e. to
- keep only zero-projection products; the
- default is to keep everything
- S -singet state density matrix (Hilbert space)
- or state vector (Liouville space)

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `exist()`, `grumble()`, `irr_sph_ten()`, `lin2lm()`, `int2str()`, `state()`, `isscalar()`, `isfield()`, `ismember()`.
