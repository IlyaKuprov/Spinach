# kernel/cce/cce_controls.m

[MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/cce/cce_controls.m)

## Purpose

Compile generalized controls, selective pulses and time-dependent exterior fields.

## Behaviour

- Optional public fields are checked by `cce_validate`; unknown, reserved, and method-inapplicable requests are rejected rather than silently ignored.

- Ideal operator events use angle one when `angle` is omitted or empty, including alongside finite controls; axis events still require an explicit finite angle.

- Coherence is raw by default on every dispatch path; explicit `normalized=true` divides by the initial coherence.
- Observable readout and error bounds both use `parameters.operators`.
- DEER species filters use the prepared model’s canonical isotope metadata.

- Expands named hahn/cpmg/xy4/xy8/udd ideal sequences and explicit piecewise-constant finite drives; amplitude and detuning are angular frequencies.
- Finite drives and ideal bath/selective rotations require generalized CCE. Simultaneous finite drives are refused rather than implicitly combined.
- Supplies piecewise segments with embedded cluster operators; exterior rotations update mean fields according to `meanfield=none|static`; static exterior ideal one-site rotations update the sampled field exactly.

- The memory budget includes retained Hamiltonians, pulses and exterior expectations for every requested time and segment, before operator allocation; scheduler accounting includes each simultaneously active worker.

- Explicit pulse-generator Hermiticity is checked relative to its norm, with no absolute unit floor; reciprocal operator/angle rescaling cannot conceal a non-Hermitian generator.

## Inputs and outputs

```matlab
segments=cce_controls(model,sites,context,schedule,parameters)
```

**Inputs**

- `model`, `sites`, `context` — physical model, included bath IDs and one global state assignment.
- `schedule` — compiled schedule, or `[]` to compile parameters.
- `parameters` — sequence/control fields, optional DEER and meanfield_radius settings; finite events specify target, placement, at, duration in seconds, axis, amplitude and detuning in rad/s.

**Outputs**

- `segments` — schedule, selected_sites, durations, hamiltonians, pulses and exterior_expectations indexed by time then segment; generalized propagation also returns raw with central densities and reducer fields.

## References

- [Source file on GitHub](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/cce/cce_controls.m)
- [cce_sequence.m](cce_sequence.md)
- [cce_cluster.m](cce_cluster.md)
