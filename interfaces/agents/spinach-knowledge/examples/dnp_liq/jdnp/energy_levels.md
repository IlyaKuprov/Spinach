# examples/dnp_liq/jdnp/energy_levels.m

- MATLAB implementation: [examples/dnp_liq/jdnp/energy_levels.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/dnp_liq/jdnp/energy_levels.m)

## What it demonstrates

A no-argument Hamiltonian illustration of how two electron-spin energy levels change from a Zeeman-dominated regime toward exchange-dominated mixing. The function deliberately exaggerates the g-factor difference to make the crossover visible; it is a spectrum/energy-diagram calculation, not a DNP polarisation simulation.

## Setup and model

Run in MATLAB with Spinach on the path. There are no function arguments or external data files. The example sets `sys.magnet=14.1` T (commented as a 600 MHz magnet), `sys.isotopes={'E','E'}`, and scalar g values `{1.9,2.1}`. The Hilbert-space basis uses `bas.formalism='zeeman-hilb'` and `bas.approximation='none'`.

Spinach builds the system and basis, then constructs the lab-frame Hamiltonian `Hz`. The exchange operator `Lj` is explicitly the isotropic dot product `Lx(1)Lx(2)+Ly(1)Ly(2)+Lz(1)Lz(2)`, with no additional prefactor in this file.

## Sweep and output

The source defines `omega_e=sys.magnet*spin('E')` and samples `omega_j` at 1000 evenly spaced points from `-3*omega_e` to `+3*omega_e`. At each point it diagonalises `Hz-omega_j*Lj`, sorts the eigenvalues, and stores the sorted energies. The plot scales both axes by `omega_e`: horizontal `omega_j/omega_e`, vertical `omega/omega_e`. Only a figure is created; the energy matrix is not returned or written to disk.

**Source-specific clarification:** the sweep diagonalises an exchange-perturbed Hamiltonian directly and sorts eigenvalues independently at each point; it does not track eigenvectors/state identity through crossings. This makes it suitable for viewing the spectrum envelope/crossing structure, not for assigning continuous state labels from line order.

**Caveats:** no relaxation, microwave irradiation, or proton is included. The deliberately large g-factor contrast is pedagogical, not an experimental parameter set.