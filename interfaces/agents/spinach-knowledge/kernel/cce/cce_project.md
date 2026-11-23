# kernel/cce/cce_project.m

## Native role and invocation

Expert conditional observed-manifold projection for compiled retained Hamiltonians. The native context requests projected generators through cce_generators.

Public simulations enter through `clexpn(spin_system,pulse_sequence,parameters,assumptions)` with a registered experiment.

Internal/expert signature: `manifolds=cce_project(spin_system,parameters)`. Engine records are not public experiment parameter structures.

[MATLAB source](../../../../../kernel/cce/cce_project.m)

## Purpose

Project an observed-spin system onto two explicit orthogonal states.

## Behaviour

- Computes manifold energies, projectors, and spin expectations in observed physical-ID tensor order.
- Optionally constructs second-order susceptibilities in the observed eigenbasis. The correction contracts ordered bath operators and vanishes only when the supplied problem makes it vanish.
- Second order is off unless enabled; only projected `cce`/`mecce` accept it, and declared gap and overlap tolerances reject unresolved projections.

## Inputs and outputs

**Inputs**

- `spin_system` — spin system returned by `cce_compile`.
- `parameters.transition.alpha`, `.beta` — normalised orthogonal state vectors. Optional `parameters.second_order.enabled`; when true, `parameters.second_order.gap_tolerance` in rad/s and `overlap_tolerance` are mandatory.

**Outputs**

- `manifolds` — alpha/beta state, projector, energy, expectations, and susceptibility records, with observed-spin system operators and projection metadata.

Validation is local to this service and its consumed numerical records. Public experiment and approximation options are validated by `cce_options` before native context dispatch.
