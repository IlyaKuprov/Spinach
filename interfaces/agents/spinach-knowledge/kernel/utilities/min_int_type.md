# kernel/utilities/min_int_type.m

- Signature: `type=min_int_type(max_val,issigned)`

## Purpose

Minimum integer data type sufficient to store the specified value. Useful in many indexing operati- ons in the Spinach kernel where double precision would be a massive overkill. Syntax: type=min_int_type(max_val,issigned)

## Physical / mathematical content

- General mathematical and infrastructure utilities. This area contains finite differences, perturbation theory, graph algorithms, spectral densities, tensor algebra, hash/report helpers, and other reusable numerical components.

## Numerical / algorithmic content

## Parameters / inputs

- max_val -maximum value that the integer
- must cover
- issigned -whether the integer needs to
- cover the negative values:
- 'signed' or 'unsigned'
- Output:
- type -Matlab data type to use

## Implementation structure

- Minimum integer data type sufficient to store the
- specified value. Useful in many indexing operati-
- ons in the Spinach kernel where double precision
- would be a massive overkill. Syntax:
- type=min_int_type(max_val,issigned)
- max_val -maximum value that the integer
- must cover
- issigned -whether the integer needs to
- cover the negative values:
- 'signed' or 'unsigned'
- Output:
- type -Matlab data type to use
