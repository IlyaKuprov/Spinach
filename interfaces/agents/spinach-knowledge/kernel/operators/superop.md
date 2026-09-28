# kernel/operators/superop.m

- Signature: `A=superop(spin_system,opspec,side)`

## Purpose

Constructs a sided product superoperator in the spherical tensor basis for left or right multiplication of a density matrix by a specified operator.

## Physical / mathematical content

- `side='left'` or `'right'` selects multiplication from that side; `'comm'` and `'acomm'` select commutation and anticommutation superoperators, respectively.
- Requires the `sphten-liouv` formalism.

## Numerical / algorithmic content

- Uses the left or right product pages of the Lie structure tables and matches source and destination basis states to assemble the result.
- Returns the superoperator in XYZ sparse format, not MATLAB's CSC format.

## Parameters / inputs

- `spin_system`: must contain basis-set information; run `basis()` before calling this function.
- `opspec`: Spinach operator specification described in Sections 2.1 and 3.3 of [this paper](http://dx.doi.org/10.1016/j.jmr.2010.11.008). It must be an integer row vector with one element per spin.
- `side`: `'left'`, `'right'`, `'comm'`, or `'acomm'`.

## Outputs

- `A`: a three-column array of row indices, column indices, and values, in that order.

## Notes

Direct calls to this general function are not usually required; use the friendlier `operator()` function. See also the [Spinach documentation](https://spindynamics.org/wiki/index.php?title=superop.m).