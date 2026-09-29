# kernel/carrier.m

- Signature: `H = carrier(spin_system,spins,operator_type)`
- Implementation: <https://github.com/IlyaKuprov/Spinach/blob/main/kernel/carrier.m>

## Contract

Builds the carrier part of the Zeeman Hamiltonian from the isotropic free-particle magnetogyric ratios and the user-specified Z-axis field, as represented in `spin_system.inter.basefrqs`. The source describes its use in rotating-frame transforms and average-Hamiltonian theory.

For each selected spin with nonzero `basefrqs`, the function adds that value times the spin's `Lz` operator, requested through `operator(spin_system,{'Lz'},{spin_index},operator_type)`. It then symmetrises the sum as `(H+H')/2` and passes it to `clean_up` with `spin_system.tols.liouv_zero`. The carrier source applies no unit conversion and does not state a unit for `basefrqs`; the coefficient is used as stored.

## Inputs and options

- `spin_system` — Spinach system object.
- `spins` — character array naming an isotope present in `spin_system.comp.isotopes` (for example, `'1H'`), or `'all'`. An isotope name selects all matching entries; `'all'` selects every spin.
- `operator_type` — optional character array. It defaults to `'comm'` and must be one of `'left'`, `'right'`, `'comm'`, or `'acomm'`. In Liouville space these request the left-product, right-product, commutation, and anticommutation superoperators, respectively.

## Output dimensions

`H` is the Hamiltonian matrix in Hilbert space or the corresponding superoperator in Liouville space. Its square dimensions follow the operator space and basis configured for `spin_system` (the Hilbert basis dimension in Hilbert space, or the Liouville basis dimension for the superoperator).

## Reference

- [Spinach Wiki: `carrier.m`](https://spindynamics.org/wiki/index.php?title=carrier.m)
