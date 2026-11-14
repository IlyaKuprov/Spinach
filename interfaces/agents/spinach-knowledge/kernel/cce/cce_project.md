# kernel/cce/cce_project.m

[MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/cce/cce_project.m)

## Purpose

Project a central register onto two explicit orthogonal states.

## Behaviour

- Optional public fields are checked by `cce_validate`; unknown, reserved, and method-inapplicable requests are rejected rather than silently ignored.

- Computes branch energies, projectors and spin expectations in central physical-ID tensor order.
- Optionally constructs second-order susceptibilities in the central eigenbasis. The correction contracts ordered bath operators and vanishes only when the supplied problem makes it vanish.
- Second order is off unless enabled; only conditional `cce`/`mecce` accept it, and declared gap and overlap tolerances reject unresolved projections.

## Inputs and outputs

```matlab
branches=cce_project(model,parameters)
```

**Inputs**

- `model` — prepared central-register model.
- `parameters.central.alpha`, `.beta` — normalized orthogonal state vectors. Optional `parameters.second_order.enabled`; when true, `parameters.second_order.gap_tolerance` in rad/s and `overlap_tolerance` are mandatory.

**Outputs**

- `branches` — alpha/beta state, projector, energy, expectations and susceptibility records, with register operators and projection metadata.

## References

- [Source file on GitHub](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/cce/cce_project.m)
- [cce_prepare.m](cce_prepare.md)
- [cce_cluster.m](cce_cluster.md)
