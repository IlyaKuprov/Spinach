# cce_state — operator states and explicit physical populations

[MATLAB source](../../../../../kernel/cce/cce_state.m)

## Native role and invocation

Expert tensor-product density preparation utility for an already constructed retained subsystem.

Public simulations enter through `clexpn(spin_system,pulse_sequence,parameters,assumptions)` with a registered experiment.

Internal/expert signature: `rho=cce_state(spin_system,states,spins)`. Engine records are not public experiment parameter structures.

Keywords: cce_state clexpn density populations trace positivity vectorisation Hilbert Liouville physical IDs

`cce_state` uses the same character descriptors and cell-row product descriptors as cce_operator for conventional unnormalised spin-density operators. Hilbert output is a sparse matrix; zeeman-liouv output is its column-major vector, not a commutator.

Explicit physical states use a cell row of local density matrices or population vectors and a matching numeric physical-ID row. Every retained spin must be specified exactly once: omitted populations are rejected, not guessed. Population vectors use descending magnetic quantum number. Matrices must be positive, Hermitian, and trace one within the declared tolerance. The output follows spin_system.cce.ids irrespective of descriptor input order, and supplied populations are not renormalised. The resolved memory budget is enforced during construction. No global full-bath state is constructed.

Operator labels and physical-ID selections, or explicit local physical densities, are validated before constructing the retained tensor product. The same positivity and trace tolerances apply to population-vector and density-matrix forms.
