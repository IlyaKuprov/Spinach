# cce_operator — canonical retained-cluster spin descriptors

[MATLAB source](../../../../../kernel/cce/cce_operator.m)

## Native role and invocation

Expert physical-spin operator constructor for an already constructed retained subsystem.

Public simulations enter through `clexpn(spin_system,pulse_sequence,parameters,assumptions)` with a registered experiment.

Internal/expert signature: `A=cce_operator(cluster_system,operators,spins)`. Engine records are not public experiment parameter structures.

Keywords: cce_operator clexpn physical IDs isotope sum product operator Zeeman commutator normalisation

`cce_operator` supports E, Lx, Ly, Lz, L+, and L- with dimensionless canonical spin normalisation. A numeric physical-ID row, isotope label, or 'all' selects a sum; matching cell rows of labels and selections specify an ordered product of sums. IDs are physical, not local indices. Unsupported labels, absent isotopes, and discarded IDs fail explicitly.

The service reuses pauli for arbitrary multiplicities and embeds in the subsystem's observed-first ordering. Hilbert output is sparse. In zeeman-liouv it returns the conventional commutator kron(I,A)-kron(A.',I) for column-major density vectors. It does not add a stock formalism or modify kernel/operator.m. The resolved model memory budget is checked before allocation. State descriptors and physical populations belong to cce_state.

All operator labels, isotope selections, and physical-ID rows are checked before retained-operator allocation, including every factor of a product-of-sums descriptor.
