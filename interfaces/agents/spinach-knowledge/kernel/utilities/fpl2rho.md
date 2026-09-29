# kernel/utilities/fpl2rho.m

## Purpose

Converts a Fokker–Planck state vector into a Liouville space state vector by averaging the spin state vector across the spatial dimensions of the sample.

## Behaviour

- Syntax: `rho=fpl2rho(rho,dims)`.
- A consistency check (`grumble`) is performed first:
  - `rho` must be a numeric array.
  - `dims` must be a row vector of real positive integers.
  - `numel(rho)` must be divisible by `prod(dims)`, i.e. it must match the space(x)spin Kronecker product; otherwise an error is raised.
- The stack size is taken as `size(rho,2)`.
- The spin dimension is exposed by reshaping the full (densified) state vector into `[size(rho,1)/prod(dims), prod(dims), stack_size]`; the comment notes there is no N-D sparse support yet.
- The spatial coordinates are averaged out via `sum(rho,2)/prod(dims)`, and the result is squeezed before being returned.

## Inputs and outputs

**Inputs**

- `rho` — Fokker–Planck state vector (numeric array).
- `dims` — spatial dimensions of the Fokker–Planck problem, a row vector of positive integers.

**Outputs**

- `rho` — Liouville space state vector.

## References

- Source: [kernel/utilities/fpl2rho.m](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/utilities/fpl2rho.m)
- Spin Dynamics Wiki: [fpl2rho.m](https://spindynamics.org/wiki/index.php?title=fpl2rho.m)
