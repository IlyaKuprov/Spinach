# kernel/pulses/bloch_siegert.m

- Signature: `[ctrl_opers,ctrl_coefs]=bloch_siegert(spin_system,ctrl_opers,ctrl_coefs)`

## Purpose

Adds the Bloch–Siegert response channel for each Cartesian control channel. The function obtains the response operators from the channel isotopes and carrier frequencies stored in `spin_system.control`, then appends each operator and the element-wise square of its matching control-amplitude vector to the input cell arrays.

## Inputs

- `spin_system` — Spinach system structure with Bloch–Siegert corrections enabled by `optimcon()`, plus its control-channel isotopes, channel indices, and carrier frequencies.
- `ctrl_opers` — cell array with one square numeric control operator per control channel.
- `ctrl_coefs` — cell array with one real vector of control coefficients (in rad/s) per channel. All coefficient vectors must have the same length.

## Outputs

- `ctrl_opers` — the input operators followed by their corresponding Bloch–Siegert response operators.
- `ctrl_coefs` — the input coefficient vectors followed by their element-wise squares, in the same channel order.

## Use and constraints

The augmented arrays are intended for the piecewise-constant propagation methods of `shaped_pulse_xy`. They are not valid for piecewise-linear interpolation: interpolating squared coefficients is not equivalent to squaring the interpolated control amplitude.

`optimcon()` keeps response operators on parallel workers; this routine rebuilds them with `bss_ops` from the channel isotopes and carrier frequencies in the returned control structure. Do not edit those settings after `optimcon()` has run, because the optimiser would then use operators built from different settings than those used by the pulse calculation.

## Reference

[Spin Dynamics Wiki: `bloch_siegert.m`](https://spindynamics.org/wiki/index.php?title=bloch_siegert.m)
