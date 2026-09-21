# kernel/operator.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/operator.m`
- Signature: `A=operator(spin_system,operators,spins,operator_type,format)`
- Total lines: 271

## Purpose

Generates Hilbert space operators or Liouville space superoperators from their human-readable descriptions. Syntax: A=operator(spin_system,operators,spins,operator_type,format)

## Physical / mathematical content

- This file belongs to the `kernel` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content

- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `numel()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Parameters / inputs

- 1. If operators is a string and spins is a string
- operators='Lz'; spins='13C';
- the function returns the sum of the corresponding single-spin operators
- (Hilbert space) or superoperators (Liouville space) on all spins of that
- type. Valid labels for states in this type of call are 'E' (identity),
- 'Lz', 'Lx', 'Ly', 'L+', 'L-', 'Tl,m' (irreducible spherical tensor, l
- and m are integers), 'CTx', 'CTy', 'CTz', 'CT+', 'CT-' (central transi-
- tion operators in the Zeeman basis). Valid labels for spins are standard
- isotope names, as well as 'electrons', 'nuclei', and 'all'.
- 2. If operators is a string and spins is a vector
- operators='Lz'; spins=[1 2 4];
- the function returns the sum of all single-spin operators (Hilbert space)
- or superoperators (Liouville space) for all spins with the specified num-
- bers. Valid labels for operators are the same as in Item 1 above.
- 3. If operators is a cell array of strings and spins is a cell array of
- numbers
- operators={'Lz','L+'}; spins={1,2};
- then a product operator (Hilbert space) or its superoperator (Liouville
- space) is produced. In the case above, Spinach will generate LzS+ in Hil-
- bert space or its specified superoperator in Liouville space. Valid la-
- bels for operators are the same as in Item 1 above.
- In Liouville space calculations, operator_type can be set to:
- 'left' -produces left side product superoperator
- 'right' -produces right side product superoperator
- 'comm' -produces commutation superoperator (default)
- 'acomm' -produces anticommutation superoperator
- In Hilbert space calculations operator_type parameter is ignored, and the
- operator itself is always returned.
- The format parameter refers to the format of the output: 'csc' returns a
- Matlab sparse matrix, 'xyz' returns a [rows, cols, vals] array.

## Outputs

- A -a CSC sparse (default) or a [rows, cols, vals] repre-
- sentation of a spin operator or superoperator.
- Notes: WARNING -a product of two commutation superoperators is NOT a com-
- mutation superoperator of a product. In Liouville space, you cannot
- generate single-spin superoperators and multiply them up.
- Note: operator caching is supported, add 'op_cache' to sys.enable array
- to enable; make sure your scratch storage is fast.

## Implementation structure

- Generates Hilbert space operators or Liouville space superoperators from
- their human-readable descriptions. Syntax:
- A=operator(spin_system,operators,spins,operator_type,format)
- 1. If operators is a string and spins is a string
- operators='Lz'; spins='13C';
- the function returns the sum of the corresponding single-spin operators
- (Hilbert space) or superoperators (Liouville space) on all spins of that
- type. Valid labels for states in this type of call are 'E' (identity),
- 'Lz', 'Lx', 'Ly', 'L+', 'L-', 'Tl,m' (irreducible spherical tensor, l
- and m are integers), 'CTx', 'CTy', 'CTz', 'CT+', 'CT-' (central transi-
- tion operators in the Zeeman basis). Valid labels for spins are standard
- isotope names, as well as 'electrons', 'nuclei', and 'all'.

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `exist()`, `grumble()`, `ismember()`, `md5_hash()`, `gcp()`, `getCurrentValueStore()`, `isKey()`, `store()`, `human2opspec()`, `superop()`, `coeffs()`, `lin2lm()`, `irr_sph_ten()`, `mults()`, `strcmp()`, `hilb2liouv()`.
