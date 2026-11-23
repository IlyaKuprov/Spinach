# cce_evaluate — one cluster of a cluster correlation expansion

`raw=cce_evaluate(spin_system,cluster,configuration,pulse_sequence,parameters)` builds the retained subsystem of one cluster with `cce_subsystem`, its generators with `cce_generators`, and calls the experiment as `pulse_sequence(spin_system,parameters,H,R,K)`. It is called by `cce_schedule`, `cce_interlace`, and `cce_partition`; public calculations enter through `clexpn`.

The returned record carries `cluster_id`, `value`, `abs_error_bound`, `valid`, and `weight`, plus the diagnostic fields that the experiment returned: generator evidence, spectral residuals, pulse timing, density, and noise moments. These records are consumed by `cce_reduce`; they are not the first output of `clexpn`.

The grumbler rejects `cce.checkpoint`, because only the scheduler owns checkpoints, and checks the retained time-history memory, including density components and cell receivers, before a subsystem is built. Projected propagation is rejected with `cce:population_transfer` when the retained Hamiltonian leaks out of either selected observed manifold by more than dimension*eps*norm(H,fro), unless the second-order treatment or the diagonal register approximation is requested. Every experiment must return its own absolute numerical estimate.

What an experiment supports (projected generators, dissipation, density output) is tested in the grumbler of that experiment.

[MATLAB source](../../../../../kernel/cce/cce_evaluate.m)
