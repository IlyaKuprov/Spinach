# kernel/cce/cce_controls.m

[MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/cce/cce_controls.m)

## Purpose

Compile generalized controls, selective pulses and time-dependent exterior fields.

## Behaviour

- Expands named hahn/cpmg/xy4/xy8/udd ideal sequences and explicit piecewise-constant finite drives; amplitude and detuning are angular frequencies.
- Finite drives and ideal bath/selective rotations require generalized CCE. Simultaneous finite drives are refused rather than implicitly combined.
- Supplies piecewise segments with embedded cluster operators; exterior rotations update mean fields according to `meanfield=none|static`; static exterior ideal one-site rotations update the sampled field exactly.

## Inputs and outputs

```matlab
segments=cce_controls(model,sites,context,schedule,parameters)
```

**Inputs**

- `model`, `sites`, `context` — physical model, included bath IDs and one global state assignment.
- `schedule` — compiled schedule, or `[]` to compile parameters.
- `parameters` — sequence/control fields, optional DEER and meanfield_radius settings; finite events specify target, placement, at, duration in seconds, axis, amplitude and detuning in rad/s.

**Outputs**

- `segments` — compiled schedule and/or time-segment Hamiltonians, pulses and exterior-field metadata.

## References

- [Source file on GitHub](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/cce/cce_controls.m)
- [cce_sequence.m](cce_sequence.md)
- [cce_cluster.m](cce_cluster.md)
