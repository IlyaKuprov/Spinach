# cce_state — operator states and explicit physical populations
Keywords: cce_state clexpn density populations trace positivity vectorisation Hilbert Liouville physical IDs

`rho=cce_state(cluster_system,states,spins)` uses the same character descriptors and cell-row product descriptors as cce_operator for conventional unnormalised spin-density operators. Hilbert output is a sparse matrix; zeeman-liouv output is its column-major vector, not a commutator.

Explicit physical states use a cell row of local density matrices or population vectors and a matching numeric physical-ID row. Every retained spin must be specified exactly once: omitted populations are rejected, not guessed. Population vectors use descending magnetic quantum number. Matrices must be positive, Hermitian and trace one within the declared tolerance. The output follows cluster_system.physical_ids irrespective of descriptor input order, and supplied populations are not renormalised. The resolved memory budget is enforced during construction. No global full-bath state is constructed and kernel/state.m is unchanged.
