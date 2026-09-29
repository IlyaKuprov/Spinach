# kernel/utilities/expdrop.m

## Purpose

Generates an exponential fall-off (drop) from a specified starting value to a specified ending value over a given duration, with a specified exponential rate and number of discretisation points.

## Behaviour

- Syntax: `drop=expdrop(from,to,duration,npoints,drop_rate)`
- The function first validates all inputs via an internal consistency-checking subfunction `grumble`, which errors out with descriptive messages if any argument fails its checks.
- The exponential parameters are computed as:
  - `B=(from-to)/(1-exp(-drop_rate*duration))`
  - `A=from-B`
- The drop is then evaluated as `drop=A+B*exp(-drop_rate*linspace(0,duration,npoints))`, returning a row vector of `npoints` values spanning the interval `[0,duration]`.
- Input validation rules enforced by `grumble`:
  - `npoints` must be a positive real integer (numeric, real, scalar, at least 1, and integer-valued).
  - `from` must be a real numeric scalar.
  - `to` must be a real numeric scalar.
  - `duration` must be a positive real numeric scalar.
  - `drop_rate` must be a positive real numeric scalar.

## Inputs and outputs

Inputs:
- `from` — the value to drop from (real scalar).
- `to` — the value to drop to (real scalar).
- `duration` — drop duration, seconds (positive real scalar).
- `npoints` — the number of discretisation points in the drop (positive integer).
- `drop_rate` — exponential drop rate, Hz (positive real scalar).

Output:
- `drop` — a row vector with the fall-off.

## References

- Source: [kernel/utilities/expdrop.m](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/utilities/expdrop.m)
- Spin Dynamics Wiki: [expdrop.m](https://spindynamics.org/wiki/index.php?title=expdrop.m)
