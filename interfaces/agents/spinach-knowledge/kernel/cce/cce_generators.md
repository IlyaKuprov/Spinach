# cce_generators — native retained generators

Syntax: [H,R,K,parameters,cache]=cce_generators(cluster_system,parameters,cache).

The subsystem is cce_subsystem output; parameters uses the native cce_options contract. Conventional generators obey L=H+1i*R+1i*K and exp(-1i*L*t), with H in rad/s and R a GKSL dissipator in s^-1; K is zero. Generalised coherent evolution supports zeeman-hilb and zeeman-liouv; conventional dissipation requires zeeman-liouv. Column-major vectorisation gives kron(I,H)-kron(H.',I). parameters.bath_jumps records have support, operator and rate; descriptors are L+, L-, Lx, Ly, Lz or explicit matrices. Only selected bath spins may be jump support. Jump support uses ordered physical IDs; a jump is included only when its whole support is retained. The raising operator L+ pumps the lower magnetic sublevel towards the upper one. This stage preserves the stage-1 bath-only support contract; observed-spin T1 jumps are not accepted by the native input schema.

Projected output is an explicitly tagged H structure: kind=cce_projected, left/right conditional bath Hamiltonians, formalism, physical_ids, dissipator and swap_matrix. R and K are empty. Pass this object directly to cce_evolution, never to stock evolution or an arbitrary pulse sequence. Swap left/right for refocusing; for projected dissipation, the swapped dissipator is swap_matrix*conj(dissipator)*swap_matrix. Projection-changing Hamiltonians/jumps and dissipative virtual-state corrections are rejected. The existing static coherent second-order susceptibility is retained.

cce_backend selects auto/direct/stock; this is an expert service construction override, not a translated legacy field. Returned parameters.cce_aux reports the retained Hamiltonian construction metadata and rounding bound. The exact generator cache is separate from the subsystem cache; changed orientation, configuration, parameters or formalism invalidates reuse. No full-bath basis is allocated, and direct construction does not repeat create/basis.

## Native dissipative, noise and ensemble contract

Conventional evolution uses `L=H+1i*R+1i*K`; GKSL rates are inverse seconds. Projected generators carry left/right branch Hamiltonians and a coherence-block dissipator, not a physical density matrix. Branch exchange transforms the dissipator with `swap_matrix*conj(dissipator)*swap_matrix'`. Ordinary and projected dissipative experiments share this service.

First-order projected Hamiltonians are conditional approximations, not a requirement that the full Hamiltonian preserve observed populations. The omitted-manifold coupling norm is recorded as `projection_leakage`. Jump operators must still preserve the selected manifolds. With nonempty GKSL jumps, projected master dynamics additionally retains the existing full-Hamiltonian manifold-invariance requirement.
