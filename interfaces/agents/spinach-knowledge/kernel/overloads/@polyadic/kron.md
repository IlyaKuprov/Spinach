# kernel/overloads/@polyadic/kron.m

- Signature: `c=kron(a,b)`

## Purpose

Kronecker product function for polyadics. Syntax: c=kron(a,b)

## Physical / mathematical content

- Polyadic tensor-product linear algebra. The emphasis is compressed operator representation, deferred algebra, and efficient Kronecker-structured manipulations.

## Numerical / algorithmic content

Core lists are extended directly only when the polyadic operand has neither prefixes nor suffixes. Otherwise, nesting preserves the complete matrix product, including multiple rectangular affixes, in either operand order. This also preserves voxel-wise flow generators when `v2fplanck` extends them into spin space.

## Parameters / inputs

- a,b -polyadic or numeric objects

## Outputs

- c -polyadic object
- This operation bundles the inputs into a nested polyadic object.

## Header notes

Both numeric and polyadic matrix factors are supported. An affixed operand is retained as a nested factor so that extension does not change its existing matrix dimensions.
