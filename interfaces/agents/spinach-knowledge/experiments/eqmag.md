# experiments/eqmag.m

- Signature: `magn=eqmag(spin_system,parameters)`

## Purpose

Compute the thermal-equilibrium molar magnetization vector at `inter.temperature` and `sys.magnet`, with the magnetic field assumed to lie along Z. Average the result over molecular orientations on the spherical grid named by `parameters.grid`.

## Calculation and output

For each grid orientation, the function rotates each spin's g-tensor, constructs the magnetic-moment operators, forms and trace-normalizes the equilibrium density matrix, and accumulates weighted expectation values. The orientation loop uses `parfor`. The return value is a real 1-by-3 vector `[Mx My Mz]` in `Na*mu_bohr`.

## Requirements and convention

- `spin_system.bas.formalism` must be `zeeman-hilb`.
- `parameters.grid` must be a nonempty character string naming a spherical averaging grid.
- Spinach uses the exchange Hamiltonian term `2*pi*J*(LxSx+LySy+LzSz)`, with `J` in Hz.
