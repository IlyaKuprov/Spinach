# cce_evolution — retained-space propagation

[MATLAB source](../../../../../kernel/cce/cce_evolution.m)

## Native role and invocation

Expert retained-space propagation kernel used by control/evaluation services. Its record contract is internal, not the clexpn public input structure.

Public simulations enter through `clexpn(spin_system,pulse_sequence,parameters,assumptions)` with a registered experiment.

Internal/expert signature: `answer=cce_evolution(cluster_system,L,coil,rho,timestep,nsteps,output)`. Engine records are not public experiment parameter structures.

L is a conventional H+1i*R+1i*K or the explicitly tagged projected structure from cce_generators. timestep is seconds, nsteps is a nonnegative integer. A coherent Hilbert state is a square density matrix; Liouville and projected bath-coherence states are column-major vectors. Non-Hermitian Hilbert generators and unsupported representation/output pairs are rejected.

Outputs: final returns the state at nsteps*timestep; trajectory includes t=0 and uses columns for Liouville/projected states and pages for Hilbert matrices; observable includes t=0 and has time-by-detector orientation; total sums states including t=0. A cell row of detectors is supported. Detection is coil(:)'*rho(:), including non-Hermitian coils. No normalisation or receiver phase is applied in evolution. Coherent propagation preserves the established eigensolver formulas; dissipative propagation uses expm(-1i*L*t). A finite pulse adds its angular-frequency control Hamiltonian before propagation. Small retained-space propagation does not rebuild a stock basis.
