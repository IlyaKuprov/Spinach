# experiments/eqmag.m

- MATLAB implementation: [experiments/eqmag.m](https://github.com/IlyaKuprov/Spinach/blob/main/experiments/eqmag.m)

- Signature: `magn=eqmag(spin_system,parameters)`

## Purpose and calculation

`eqmag` returns the thermal-equilibrium molar magnetisation vector averaged over the molecular orientations in a spherical grid. The function assumes the magnetic field specified by `spin_system.inter.magnet` is along Z and uses the temperature in `spin_system.rlx.temperature`. It requires one substance and `spin_system.bas.formalism='zeeman-hilb'`; segmented systems raise `Spinach:eqmag:multipleSubstances` because mixture molar normalisation is not defined by this interface; `parameters.grid` must be a nonempty character string naming a grid available to the Spinach grid loader.

The caller does not prepare magnetisation operators: the routine obtains each spin's `Lx`, `Ly`, and `Lz` operators and g-tensor, rotates the g-tensors for each grid orientation, and evaluates the corresponding equilibrium density matrix. The grid weights average the orientation-dependent traces. This is an equilibrium calculation, not a driven or time-resolved pulse sequence.

## Output and convention

`magn` is a real 1-by-3 vector `[Mx My Mz]` in `Na*mu_bohr` units (Avogadro's constant times the Bohr magneton). The source also records Spinach's exchange-coupling convention: `2*pi*J*(LxSx+LySy+LzSz)`, with `J` in Hz.

[Source page](https://spindynamics.org/wiki/index.php?title=eqmag.m)
