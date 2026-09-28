# kernel/overloads/@polyadic/validate.m

- Signature: `validate(p)`

## Purpose

Checks the structure and factor dimensions of a polyadic object, raising an error when a representation invariant is not met.

## Physical / mathematical content

This is a representation validator, not a matrix operation: it checks that core terms, prefix factors, and suffix factors can form compatible matrix products.

## Numerical / algorithmic content

Dimension totals for each core term are obtained by multiplying the row and column sizes of its stored core matrices. The routine also reports a warning when the number of core terms or boundary factors exceeds 100.

## Parameters / inputs

- p -a polyadic object

## Implementation structure

- Requires p to be a polyadic and its top-level cores, prefix, and suffix fields to be cell arrays.
- Checks that each core term is a cell array and that its entries and boundary factors are numeric; the implementation also calls validate recursively for entries identified as polyadic.
- Checks that all buffered core terms have the same total row and column dimensions.
- Checks that prefix-to-core, core-to-suffix, and successive prefix/suffix dimensions match.
- Warns if the number of core terms, prefix factors, or suffix factors is greater than 100.
