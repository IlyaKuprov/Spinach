# kernel/utilities/spden.m

- Signature: `J=spden(L,D,omega)`

## Purpose

Computes the Lorentzian spectral density for rotational diffusion at the specified frequency.

## Physical / mathematical content

The correlation time and spectral density are

- `tau_c=1/(L*(L+1)*D)`
- `J=(tau_c/(2*L+1))/(1+(tau_c*omega)^2)`

## Numerical / algorithmic content

The function checks the inputs, calculates `tau_c`, then calculates `J`.

## Parameters / inputs

- `L` — spherical rank; a positive real integer. Use 2 for common NMR mechanisms such as dipolar relaxation.
- `D` — rotational diffusion coefficient in s⁻¹; a positive real scalar.
- `omega` — frequency in rad/s; a real scalar.

## Outputs

- `J` — spectral density function value.

## Implementation structure

Input consistency is enforced by the local `grumble` function.

Source: <https://spindynamics.org/wiki/index.php?title=spden.m>

ilya.kuprov@weizmann.ac.il