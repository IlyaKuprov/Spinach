# kernel/utilities/oscillator.m

## Purpose

Builds 1D harmonic oscillator infrastructure: the oscillator Hamiltonian, the position operator, and the coordinate grid, given force constant, particle mass, gravitational acceleration, and discretisation settings.

## Behaviour

- Syntax: `[H_oscl,X_oscl,xgrid]=oscillator(parameters)`.
- Validates all inputs via an internal `grumble` function, which errors if any field is missing or invalid.
- Builds the second derivative operator as `((n_points-1)/box_size)^2*fdmat(n_points,5,2)`, i.e. a 5-point finite-difference second derivative scaled to the physical grid.
- Constructs `xgrid` with `linspace(-box_size/2, box_size/2, n_points)` as a column vector.
- Builds `X_oscl` as a sparse diagonal matrix with `xgrid` on the main diagonal (`spdiags(xgrid,0,n_points,n_points)`).
- Applies zero boundary conditions by zeroing the first and last rows and columns of the second derivative operator.
- Assembles the Hamiltonian as `-(1/(2*par_mass))*d2_dx2 + (frc_cnst/2)*X_oscl^2 + par_mass*grv_cnst*X_oscl`.
- Gravitation is directed along the X axis.
- The reduced Planck constant is set to 1 J*s in the kinetic energy term, so the Hamiltonian is also the generator of time evolution in rad/s, as in `expm(-1i*H_oscl*t)`.

## Inputs and outputs

Inputs (fields of `parameters`):

- `parameters.frc_cnst` — force constant, N/m; must be a positive real scalar.
- `parameters.par_mass` — particle mass, kg; must be a positive real scalar.
- `parameters.grv_cnst` — gravitational acceleration, m/s^2; must be a real scalar (may be negative or zero).
- `parameters.n_points` — number of discretisation points; must be a positive real integer.
- `parameters.box_size` — oscillator box size, m; must be a positive real scalar.

Outputs:

- `H_oscl` — oscillator Hamiltonian, Joules with hbar=1.
- `X_oscl` — oscillator X operator, m.
- `xgrid` — X coordinate grid, m.

## References

- Source: [kernel/utilities/oscillator.m](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/utilities/oscillator.m)
- Spinach Wiki: [oscillator.m](https://spindynamics.org/wiki/index.php?title=oscillator.m)
