# kernel/utilities/min_int_type.m

- Signature: `type=min_int_type(max_val,issigned)`

## Purpose

Minimum integer data type sufficient to store the specified value. Useful in many indexing operations in the Spinach kernel where double precision would be a massive overkill. Syntax: type=min_int_type(max_val,issigned)

## Physical / mathematical content

- Selects the smallest built-in MATLAB signed or unsigned integer class whose positive range covers `max_val`.

## Numerical / algorithmic content

- Compares `max_val` with `intmax` for the signed or unsigned 8-, 16-, 32-, and 64-bit classes, returning the first fit or raising an error if none fits.

## Parameters / inputs

- max_val -maximum value that the integer
- must cover
- issigned -whether the integer needs to
- cover the negative values:
- 'signed' or 'unsigned'
- Output:
- type -Matlab data type to use

## Implementation structure

- Requires a positive real integer `max_val` and `issigned` equal to `signed` or `unsigned`, then checks the corresponding integer classes in increasing width.
