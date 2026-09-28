# kernel/optimcon/distortions/no_dist.m

- Signature: `[w,J]=no_dist(w)`

## Purpose

Returns the waveform unchanged. If a second output is requested, it returns the unit Jacobian for the vectorised waveform.

## Parameters / inputs

- `w` — a numerical waveform array. The function requires its elements to be real.

## Outputs

- `w` — the input waveform, unchanged.
- `J` — if requested, a sparse identity matrix of size `numel(w)` by `numel(w)`.

## Implementation

The function checks that `w` is numeric and real, raising an error otherwise. It constructs `J` with `speye(numel(w))` only when a second output is requested.

## Reference

- [Spinach documentation for `no_dist.m`](https://spindynamics.org/wiki/index.php?title=no_dist.m)