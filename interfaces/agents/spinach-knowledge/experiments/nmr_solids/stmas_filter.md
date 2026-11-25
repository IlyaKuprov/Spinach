# STMAS transition support projector

`P=stmas_filter(spin_system,spin_idx,level_pairs)` selects explicit Zeeman
density-matrix elements of one spin, rather than an entire coherence order.
Use a complete single-substance `zeeman-liouv` basis. Each row of
`level_pairs` is an ordered density-matrix (row, column) pair in descending
magnetic-quantum-number order. For spin 3/2, `[1 2; 3 4]` selects both
positive-order satellites and `[2 3]` selects the positive-order central
transition. Swap columns to select the conjugate coherences.

Multiplying a column-vectorised spin state by `P` preserves each selected
element unweighted, including off-diagonal spectator indices. The filter
is a projector, not a raising operator: it applies no ladder coefficients
and does not mix coherences. For a rotor-space state, lift this spin-space
projector with the identity on the spatial factor.
