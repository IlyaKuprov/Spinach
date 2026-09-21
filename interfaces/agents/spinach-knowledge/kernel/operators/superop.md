# kernel/operators/superop.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/operators/superop.m`
- Signature: `A=superop(spin_system,opspec,side)`
- Total lines: 212

## Purpose

Sided product superoperator in the spherical tensor basis set. Returns superoperators corresponding to right or left multiplication of a den- sity matrix by a user-specified operator. Syntax: A=superop(spin_system,opspec,side)

## Physical / mathematical content

- Operator-construction utilities. They build bases and irreducible tensor representations for spin, bosonic, and transition operators.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Parameters / inputs

- opspec -Spinach operator specification described in Sections 2.1
- and 3.3 of the following paper:
- side -'left' or 'right' causes the function to return a product
- superoperator corresponding to a product from that side;
- 'comm' or 'acomm' results in commutation and anticommuta-
- tion superoperator respectively.

## Outputs

- A -a three-column array of row indices (first column),
- column indices (second column) and values (third column).
- Note: this is a very general function to which direct calls are not
- usually required -please use the (much friendlier) operator()
- function.
- Note: the superoperator is returned in XYZ sparse format, which is
- different from Matlab's CSC format.

## Implementation structure

- Sided product superoperator in the spherical tensor basis set. Returns
- superoperators corresponding to right or left multiplication of a den-
- sity matrix by a user-specified operator. Syntax:
- A=superop(spin_system,opspec,side)
- opspec -Spinach operator specification described in Sections 2.1
- and 3.3 of the following paper:
- side -'left' or 'right' causes the function to return a product
- superoperator corresponding to a product from that side;
- 'comm' or 'acomm' results in commutation and anticommuta-
- tion superoperator respectively.
- A -a three-column array of row indices (first column),
- column indices (second column) and values (third column).

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `strcmp()`, `grumble()`, `unit_oper()`, `active_spins()`, `opspec()`, `squeeze()`, `ismember()`, `from()`, `coeff()`, `true()`, `and()`, `basis_cols()`, `source_subsp()`, `destin_subsp()`, `isequal()`, `source_subsp_idx()`.
