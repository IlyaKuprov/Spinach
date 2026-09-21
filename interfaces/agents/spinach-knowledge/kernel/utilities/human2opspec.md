# kernel/utilities/human2opspec.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/utilities/human2opspec.m`
- Signature: `[opspecs,coeffs]=human2opspec(spin_system,operators,spins)`
- Total lines: 374

## Purpose

Converts user-friendly descriptions of spin states and operators into the formal description (opspec) used by Spinach kernel. Syntax: [opspecs,coeffs]=human2opspec(spin_system,operators,spins)

## Physical / mathematical content

- General mathematical and infrastructure utilities. This area contains finite differences, perturbation theory, graph algorithms, spectral densities, tensor algebra, hash/report helpers, and other reusable numerical components.
- The effective hardware model is a weakly anharmonic oscillator. Duffing nonlinearity breaks equal level spacing and allows qubit-like addressability within a truncated bosonic ladder.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `numel()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Parameters / inputs

- 1. If both inputs are strings
- [opspecs,coeffs]=human2opspec(spin_system,'Lz','13C')
- the function returns a list of single-spin opspecs for all spins with the
- specified name. In the example above, the list of Lz operator specificati-
- ons for all 13C nuclei in the system would be returned. Valid labels for
- states in this type of call are 'E' (identity), 'Lz', 'Lx', 'Ly', 'L+',
- 'L-', 'Tl,m' (irreducible spherical tensor, l and m are integers), 'CTx',
- 'CTy', 'CTz', 'CT+','CT-' (central transition operators in the Zeeman ba-
- sis). Valid labels for spins are standard isotope names, as well as 'elec-
- trons', 'nuclei', and 'all'.
- 2. If one input is a string and the other is a vector
- [opspecs,coeffs]=human2opspec(spin_system,'Lz',[1 2 4])
- the function returns a list of single-spin opspecs for all spins with the
- specified number. In the example above, the list of Lz operator specifica-
- tions for all 13C nuclei in the system would be returned. Valid labels for
- operators are the same as in Item 1 above.
- 3. If the two inputs are a cell array of strings and a cell array of num-
- bers, a product operator specification is produced
- [opspecs,coeffs]=human2opspec(spin_system,{'Lz','Ly'},{1,2})
- would return the Lz(x)Ly product operator specification with Lz on spin 1
- and Ly on spin 2. Valid labels for operators are the same as in Item 1.

## Outputs

- opspecs -Spinach operator specification: a cell array of
- row vectors specifying which operator enters the
- Kronecker product for which spin.
- coeffs -coefficient with which each of the Kronecker pro-
- ducts enters the overall sum.
- Notes: direct calls to this function are not necessary, use operator.m and
- state.m functions instead.

## Implementation structure

- Converts user-friendly descriptions of spin states and operators into the
- formal description (opspec) used by Spinach kernel. Syntax:
- [opspecs,coeffs]=human2opspec(spin_system,operators,spins)
- 1. If both inputs are strings
- [opspecs,coeffs]=human2opspec(spin_system,'Lz','13C')
- the function returns a list of single-spin opspecs for all spins with the
- specified name. In the example above, the list of Lz operator specificati-
- ons for all 13C nuclei in the system would be returned. Valid labels for
- states in this type of call are 'E' (identity), 'Lz', 'Lx', 'Ly', 'L+',
- 'L-', 'Tl,m' (irreducible spherical tensor, l and m are integers), 'CTx',
- 'CTy', 'CTz', 'CT+','CT-' (central transition operators in the Zeeman ba-
- sis). Valid labels for spins are standard isotope names, as well as 'elec-

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `ischar()`, `cellfun()`, `strncmp()`, `strcmp()`, `spin_numbers()`, `vertcat()`, `cell2mat()`, `spins()`, `iscell()`, `opspecs()`, `opspecs_a()`, `opspecs_b()`, `ct2ist()`, `regexp()`, `textscan()`.
