# kernel/operators/superop.m

- Signature: `A=superop(spin_system,opspec,side)`

## Purpose

Sided product superoperator in the spherical tensor basis set. Returns superoperators corresponding to right or left multiplication of a den- sity matrix by a user-specified operator. Syntax: A=superop(spin_system,opspec,side)

## Physical / mathematical content

- Operator-construction utilities. They build bases and irreducible tensor representations for spin, bosonic, and transition operators.

## Numerical / algorithmic content

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
