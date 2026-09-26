# kernel/state.m

- Signature: `rho=state(spin_system,states,spins,method)`

## Purpose

Generates Hilbert space density matrices and Liouville space state vectors from their human-readable descriptions. Syntax: rho=state(spin_system,states,spins,method)

## Physical / mathematical content

## Numerical / algorithmic content

- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.

## Parameters / inputs

- 1. If states is a string and spins is a string
- states='Lz'; spins='13C';
- the function returns the sum of the corresponding single-spin densi-
- ty matrices (Hilbert space) or state vectors (Liouville space) on
- all spins of that type. Valid labels for states in this type of call
- are 'E' (identity), 'Lz', 'Lx', 'Ly', 'L+', 'L-', 'Tl,m' (irreduci-
- ble spherical tensor, l and m are integers), 'CTx', 'CTy', 'CTz',
- 'CT+','CT-' (central transition operators in the Zeeman basis). Va-
- lid labels for spins are standard isotope names, as well as 'elect-
- rons', 'nuclei', and 'all'.
- 2. If states is a string and spins is a vector
- states='Lz'; spins=[1 2 4];
- the function returns the sum of all single-spin density matrices
- (Hilbert space) or state vectors (Liouville space) for all spins
- with the specified numbers. Valid labels for states are the same as
- in Item 1 above.
- 3. If states is a cell array of strings and spins is a cell array
- of numbers:
- states={'Lz','L+'}; spins={1,2};
- then a product state density matrix (Hilbert space) or state vector
- (Liouville space) is produced. In the case above, Spinach will gene-
- rate LzS+ density matrix in Hilbert space or its state vector in Li-
- ouville space. Valid labels for operators are the same as in Item 1
- above.
- 4. For wavefunction formalism, states must be specified as an array
- of projection quantum numbers on all spins; in that case only two
- arguments are needed, for example, in a {'1H','1H','14N'} system:
- psi=state(spin_system,[-1/2 1/2 0])
- Method argument has the following effect in sphten-liouv formalism:
- 'cheap' -the state vector is generated without
- normalisation. For very large spin sys-
- tems this is much faster
- 'exact' -exact state vector with correct normalisation,
- this is the default when the last argument is
- skipped in the function call
- 'chem' -the exact state vector weighted with the
- concentrations specified in inter.chem.concs
- field under chemical kinetics parameters
- This option is ignored in zeeman-hilb and zeeman-liouv formalisms
- because there are no cheap shortcuts and kinetics is not available.

## Outputs

- rho -a Hilbert space density matrix or a Liouville
- space state vector

## Implementation structure

- Generates Hilbert space density matrices and Liouville space state
- vectors from their human-readable descriptions. Syntax:
- rho=state(spin_system,states,spins,method)
- 1. If states is a string and spins is a string
- states='Lz'; spins='13C';
- the function returns the sum of the corresponding single-spin densi-
- ty matrices (Hilbert space) or state vectors (Liouville space) on
- all spins of that type. Valid labels for states in this type of call
- are 'E' (identity), 'Lz', 'Lx', 'Ly', 'L+', 'L-', 'Tl,m' (irreduci-
- ble spherical tensor, l and m are integers), 'CTx', 'CTy', 'CTz',
- 'CT+','CT-' (central transition operators in the Zeeman basis). Va-
- lid labels for spins are standard isotope names, as well as 'elect-
