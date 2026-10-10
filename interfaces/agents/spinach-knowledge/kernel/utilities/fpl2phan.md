# kernel/utilities/fpl2phan.m

## Purpose

Returns the image painted within the Fokker-Planck vector by a user-specified spin state, effectively projecting a Fokker-Planck state vector onto an observable to obtain a spatial phantom.

Source: <https://github.com/IlyaKuprov/Spinach/blob/main/kernel/utilities/fpl2phan.m>

## Behaviour

- Calls the internal consistency checker `grumble(rho,coil,dims)` before any computation.
- Reshapes `rho` into a matrix of size `[numel(coil) prod(dims)]`, exposing the spin dimension as the first axis.
- Computes the observable as `phan=coil'*rho`, i.e. the conjugate transpose of `coil` applied to the reshaped state.
- Reshapes the result back to `dims`, producing a spatial array.

Consistency enforcement (`grumble`):

- `rho` must be numeric and a column vector, otherwise errors with `'rho must be a column vector.'`.
- `coil` must be numeric and a column vector, otherwise errors with `'coil must be a column vector.'`.
- `dims` must be numeric, a row vector, real, with all elements at least 1 and integer-valued, otherwise errors with `'dims must be a row vector of real positive integers.'`.
- `numel(rho)` must equal `prod(dims)*numel(coil)`, otherwise errors with `'numel(rho) does not match space(x)spin Kronecker product.'`.

## Inputs and outputs

Syntax:

```matlab
phan=fpl2phan(rho,coil,dims)
```

Inputs:

- `rho` — state vector in Fokker-Planck space; numeric column vector with `numel(rho)` equal to `prod(dims)*numel(coil)`.
- `coil` — observable state vector in Liouville space; numeric column vector.
- `dims` — spatial dimensions of the Fokker-Planck problem; a row vector of real positive integers.

Output:

- `phan` — the image painted by the specified state; an array shaped according to `dims`.

## References

- Spinach Wiki: <https://spindynamics.org/wiki/index.php?title=fpl2phan.m>
- Source file: <https://github.com/IlyaKuprov/Spinach/blob/main/kernel/utilities/fpl2phan.m>
