# kernel/utilities/oscillator.m

- Signature: `[H_oscl,X_oscl,xgrid]=oscillator(parameters)`

## Purpose

Construct a one-dimensional harmonic oscillator Hamiltonian, coordinate operator, and coordinate grid. Gravitation acts along the X axis; finite-difference derivative operators are used.

## Parameters

- `parameters.frc_cnst` — force constant, N/m; a positive real scalar.
- `parameters.par_mass` — particle mass, kg; a positive real scalar.
- `parameters.grv_cnst` — gravitational acceleration, m/s^2; a real scalar.
- `parameters.n_points` — number of discretization points; a positive integer scalar.
- `parameters.box_size` — oscillator box size, m; a positive real scalar.

All five fields are required and checked by `grumble`.

## Outputs

- `H_oscl` — oscillator Hamiltonian, Joules with `hbar=1`. The reduced Planck constant is set to 1 J*s in the kinetic energy term, so the Hamiltonian is also the generator of time evolution in rad/s, as in `expm(-1i*H_oscl*t)`.
- `X_oscl` — oscillator X operator, m.
- `xgrid` — X coordinate grid, m.

## Construction

The coordinate grid is `linspace(-parameters.box_size/2,parameters.box_size/2,parameters.n_points)'`, and `X_oscl` is its diagonal operator, constructed with `spdiags`. The second-derivative operator is `((parameters.n_points-1)/parameters.box_size)^2*fdmat(parameters.n_points,5,2)`; its first and last rows and columns are set to zero for boundary handling. The Hamiltonian combines kinetic, harmonic, and gravitational terms:

`H_oscl=-(1/(2*parameters.par_mass))*d2_dx2+(parameters.frc_cnst/2)*X_oscl^2+parameters.par_mass*parameters.grv_cnst*X_oscl`.

Source contact: ilya.kuprov@weizmann.ac.il. [Function reference](https://spindynamics.org/wiki/index.php?title=oscillator.m).