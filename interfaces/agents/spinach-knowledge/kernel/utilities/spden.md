# kernel/utilities/spden.m

## Purpose

`spden` returns the value of the Lorentzian spectral density function for rotational diffusion at a user-specified frequency. It is used in relaxation-theory calculations involving rotational diffusion.

## Behaviour

The function is called as `J=spden(L,D,omega)`. It first runs an internal consistency check (`grumble`) on the inputs, then computes the rotational correlation time

`tau_c = 1/(L*(L+1)*D)`

and evaluates the spectral density as

`J = (tau_c/(2*L+1))/(1+(tau_c*omega)^2)`.

Input validation enforces:

- `L` must be a positive real integer (numeric, scalar, real, `L >= 1`, integer-valued); otherwise the error `L must be a positive real integer.` is raised.
- `D` must be a positive real scalar; otherwise the error `D must be a positive real scalar.` is raised.
- `omega` must be a real scalar; otherwise the error `omega must be a real scalar.` is raised.

The header comment notes that `L = 2` should be used for common NMR mechanisms such as dipolar relaxation.

## Inputs and outputs

Inputs:

- `L` — spherical rank; use 2 for common NMR mechanisms such as dipolar relaxation.
- `D` — rotational diffusion coefficient, s^-1.
- `omega` — frequency, rad/s.

Output:

- `J` — spectral density function value.

## References

- Source: [Spinach repository — kernel/utilities/spden.m](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/utilities/spden.m)
- Spinach Wiki: [spden.m](https://spindynamics.org/wiki/index.php?title=spden.m)
