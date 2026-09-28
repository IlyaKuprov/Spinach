# kernel/operator.m

- Signature: `A=operator(spin_system,operators,spins,operator_type,format)`

## Purpose

Construct an operator or superoperator from the requested single-spin operators and spins, in the active Spinach formalism.

## Physical / mathematical content

A string `operators` with string `spins` requests a sum over spins of that type; paired operator and spin cells specify a product. In Hilbert space the function returns the operator itself. In Liouville space, `operator_type='comm'` requests the commutation superoperator and `operator_type='acomm'` the anticommutation superoperator. The Hilbert-space branch ignores `operator_type`.

## Numerical / algorithmic content

The input specification is parsed by `human2opspec`; operator terms are constructed and combined as a sum or product. In spherical-tensor Liouville formalism, coefficients weight the superoperator terms. The Zeeman-formalism terms are assembled using Kronecker products, and term construction is parallelized.

The default output format is `csc`, a complex sparse square matrix: its dimension is the number of basis rows for `sphten-liouv`, `prod(mults)` for `zeeman-wavef` and `zeeman-hilb`, and `prod(mults)^2` for `zeeman-liouv`. Format `xyz` returns a three-column `[rows,cols,vals]` array. Operator caching is used only when `op_cache` is enabled and a worker `ValueStore` is available; without a parallel pool it is skipped.

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

The source validates the basis, permitted operator/spin input forms, cell lengths and contents, spin indices, `operator_type` type, and the `format` value (`csc` or `xyz`). It requires a character `operator_type` but does not check it against the documented `comm` and `acomm` choices. Cache access requires `op_cache` in `spin_system.enable` and an available parallel-worker `ValueStore`; otherwise the routine computes without using the cache.
