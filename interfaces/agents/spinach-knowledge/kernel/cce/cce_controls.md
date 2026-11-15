# kernel/cce/cce_controls.m

[MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/cce/cce_controls.m)

## Purpose

Compile generalized controls, selective pulses and time-dependent exterior fields.

## Behaviour

- Pulse-position and ordering arrays are checked against `numerics.max_memory_bytes`, when supplied, before allocation; driven propagation also includes retained time-dependent output records in its existing segment-storage preflight.

DEER frequency windows require axial canonical lab-frame Zeeman rows at every candidate site after species/site filtering. A transverse row raises `cce_controls:deer_axial` with the physical site ID; it is never silently replaced by its longitudinal coefficient. Existing signed axial frequencies, strict boundaries, and explicit empty-selection behavior are unchanged.

A finite `meanfield_radius` requires finite real observed and exterior coordinates; missing geometry raises `cce:coordinates` before distance filtering. The radius is nonnegative, and `Inf` retains the unbounded default without requiring coordinates. Included cluster sites need no coordinates for this cutoff, and an empty exterior requires no distance calculation.

- Optional public fields are checked by `cce_validate`; unknown, reserved, and method-inapplicable requests are rejected rather than silently ignored.

- Ideal operator events use angle one when `angle` is omitted or empty, including alongside finite controls; axis events still require an explicit finite angle.

- Signal is raw by default on every dispatch path; explicit `normalized=true` divides by the initial signal.
- Signal detection and error bounds both use `parameters.coil`.
- DEER species filters use the prepared model’s canonical isotope metadata.

- Expands named hahn/cpmg/xy4/xy8/udd ideal sequences and explicit piecewise-constant finite drives; amplitude and detuning are angular frequencies.
- Finite drives and ideal bath/selective rotations require generalized CCE. Simultaneous finite drives are refused rather than implicitly combined.
- Supplies piecewise segments with embedded cluster operators; exterior rotations update mean fields according to `meanfield=none|static`; static exterior ideal one-site rotations update the sampled field exactly.

- The memory budget includes retained Hamiltonians, pulses and exterior expectations for every requested time and segment, before operator allocation; scheduler accounting includes each simultaneously active worker.

- Explicit pulse-generator Hermiticity is checked relative to its norm, with no absolute unit floor; reciprocal operator/angle rescaling cannot conceal a non-Hermitian generator.

## Inputs and outputs

```matlab
segments=cce_controls(model,sites,spin_config,schedule,parameters)
```

**Inputs**

- `model`, `sites`, `spin_config` — physical model, included bath IDs and one global state assignment.
- `schedule` — compiled schedule, or `[]` to compile parameters.
- `parameters` — sequence/control fields, optional DEER and meanfield_radius settings; finite events specify target, placement, at, duration in seconds, axis, amplitude and detuning in rad/s.

**Outputs**

- `segments` — schedule, selected_sites, durations, hamiltonians, pulses and exterior_expectations indexed by time then segment; generalized propagation also returns raw with observed densities and reducer fields.

A single-site ideal operator pulse wholly outside the retained cluster is skipped when `meanfield=none`, including driven observed-spin-only baselines. Exterior operator pulses remain rejected with static mean fields, and partially retained/multisite operator targets are not silently discarded.

## References

- [Source file on GitHub](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/cce/cce_controls.m)
- [cce_sequence.m](cce_sequence.md)
- [cce_cluster.m](cce_cluster.md)

Finite-control gCCE uses `spin_config.rho0` both for the propagated density and for normalized observable/signal references, including the observable error bound. Parameter-level `rho0` does not override a supplied spin_config.
