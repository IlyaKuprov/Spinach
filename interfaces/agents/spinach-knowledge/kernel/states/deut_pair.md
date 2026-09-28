# kernel/states/deut_pair.m

- Signature: `[S,T,Q,Tc,Qc]=deut_pair(spin_system,spin_a,spin_b,options)`

## Purpose

All possible states of a spin-1 pair, classified by the total spin into singlet, triplet, and quintet. Syntax: [S,T,Q,Tc,Qc]=deut_pair(spin_system,spin_a,spin_b,options)

## Physical / mathematical content

- For two spin-1 particles, the combined total-spin manifolds are singlet (S), triplet (T), and quintet (Q).

## Parameters / inputs

- spin_a -the number of the first spin
- spin_b -the number of the second spin
- options.dephasing -set to 1 to eliminate the states that are
- not stationary under the Zeeman Hamilto-
- nian of two inequivalent spins, i.e. to
- keep only zero-projection products; the
- default is to keep everything

## Outputs

- S -singlet state density matrix (Hilbert space)
- or state vector (Liouville space)
- T -triplet state density matrices (Hilbert space)
- or state vectors (Liouville space), ordered in
- a cell array as {T+,T0,T-}
- Q -quintet state density matrices (Hilbert space)
- or state vectors (Liouville space), ordered in
- a cell array as {Q++,Q+,Q0,Q-,Q--}
- Tc -coherences between triplet states:
- {T0 -> T-, T+ -> T0, T- -> T0, T0 -> T+}
- Qc -coherences between quintet states:
- {Q- -> Q--, Q0 -> Q-, Q+ -> Q0, Q++ -> Q+, ...
-  Q-- -> Q-, Q- -> Q0, Q0 -> Q+, Q+ -> Q++}
- The coherence labels above are Zeeman-state products, not irreducible spherical tensors.
