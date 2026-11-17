# cce_generators — native retained generators

## Native role and invocation

Internal generator construction for one retained subsystem, called by the native context evaluator.

Public simulations enter through `clexpn(spin_system,pulse_sequence,parameters,assumptions)` with a registered experiment.

Internal/expert signature: `[H,R,K,parameters,cache]=cce_generators(cluster_system,parameters,cache)`. Engine records are not public experiment parameter structures.

Syntax: [H,R,K,parameters,cache]=cce_generators(cluster_system,parameters,cache).

The subsystem is cce_subsystem output; parameters uses the native cce_options contract. Conventional generators obey L=H+1i*R+1i*K and exp(-1i*L*t), with H in rad/s and R a GKSL dissipator in s^-1; K is zero. Generalised coherent evolution supports zeeman-hilb and zeeman-liouv; conventional dissipation requires zeeman-liouv. Column-major vectorisation gives kron(I,H)-kron(H.',I). parameters.bath_jumps records have support, operator and rate; descriptors are L+, L-, Lx, Ly, Lz or explicit matrices. Jump support uses ordered observed or selected bath IDs; an active jump is included only when its whole support is retained. The raising operator L+ pumps the lower magnetic sublevel towards the upper one. Generalised Liouville evolution accepts observed-spin population-transfer jumps; projected evolution requires invariant selected manifolds.

Projected output is an explicitly tagged H structure: kind=cce_projected, left/right conditional bath Hamiltonians, formalism, physical_ids, dissipator and swap_matrix. R and K are empty. Pass this object directly to cce_evolution, never to stock evolution or an arbitrary pulse sequence. Swap left/right for refocusing; for projected dissipation, the swapped dissipator is swap_matrix*conj(dissipator)*swap_matrix. Projection-changing Hamiltonians/jumps and dissipative virtual-state corrections are rejected. The existing static coherent second-order susceptibility is retained.

cce_backend selects auto/direct/stock; this is an expert service construction override. Returned parameters.cce_aux reports the retained Hamiltonian construction metadata and rounding bound. The exact generator cache is separate from the subsystem cache; changed orientation, configuration, parameters or formalism invalidates reuse. No full-bath basis is allocated, and direct construction does not repeat create/basis.

## Native dissipative, noise and ensemble contract

Conventional evolution uses `L=H+1i*R+1i*K`; GKSL rates are inverse seconds. Projected generators carry left/right branch Hamiltonians and a coherence-block dissipator, not a physical density matrix. Branch exchange transforms the dissipator with `swap_matrix*conj(dissipator)*swap_matrix'`. Ordinary and projected dissipative experiments share this service.

First-order projected Hamiltonians are conditional approximations, not a requirement that the full Hamiltonian preserve observed populations. The omitted-manifold coupling norm is recorded as `projection_leakage`. Jump operators must still preserve the selected manifolds. With nonempty GKSL jumps, projected master dynamics additionally retains the existing full-Hamiltonian manifold-invariance requirement.

Coherent projected generators additionally expose the compiled Hamiltonian round-off estimate and full retained Hamiltonian Frobenius norm for experiment-owned propagation error accounting.

Zero-rate records are omitted before support and manifold-closure checks. An active jump split across retained and exterior spins is rejected with cce:jumpSupport. A wholly exterior active jump is omitted only when it does not drive a declared static exterior field; otherwise cce:exteriorJump is raised. The projected jump and loss-operator leakage checks retain their relative operator-norm bounds, without an absolute scale floor.

Resident jump and superoperator storage is checked before Hamiltonian construction and projected-manifold closure. An insufficient allocation budget is therefore rejected before building matrices or examining leakage, on both projected and generalized direct-service paths. The existing storage formula and physical jump algebra are unchanged.
