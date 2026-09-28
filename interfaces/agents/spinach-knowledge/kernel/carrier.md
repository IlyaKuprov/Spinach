# kernel/carrier.m

- Signature: `H=carrier(spin_system,spins,operator_type)`

## Purpose

Returns the "carrier" Hamiltonian—the part of the Zeeman interaction Hamiltonian corresponding to particles at the Zeeman frequency prescribed by their isotropic free-particle magnetogyric ratio and the user-specified Z-axis magnetic field. This Hamiltonian is used in rotating-frame transforms and average Hamiltonian theories.

## Physical / mathematical content

The carrier Hamiltonian is a sum of the selected spins' base frequencies multiplied by their `Lz` operators.

## Numerical / algorithmic content

- Selects all spins when `spins` is `'all'`; otherwise, selects spins whose isotope matches `spins`.
- Adds a frequency-weighted `Lz` operator for each selected spin with a nonzero base frequency.
- Symmetrizes the result as `(H+H')/2` and applies `clean_up` with `spin_system.tols.liouv_zero`.

## Parameters / inputs

- `spin_system` - the spin system.
- `spins` - a string specifying the isotope, e.g. `'1H'`; use `'all'` to select all spins.
- `operator_type` - in Liouville space, `'left'` produces a left-side product superoperator, `'right'` a right-side product superoperator, `'comm'` a commutation superoperator (default), and `'acomm'` an anticommutation superoperator. In Hilbert space this parameter is ignored.

## Outputs

- `H` - a Hamiltonian (Hilbert space) or its superoperator of the specified type (Liouville space).

## Implementation structure

The function validates the inputs, preallocates `H`, selects the spins, accumulates their carrier terms, and cleans up the symmetrized result.

<https://spindynamics.org/wiki/index.php?title=carrier.m>
