# kernel/overloads/@ttclass/numel.m

- Signature: `n=numel(tt)`

## Purpose

Returns the number of elements in the matrix represented by a tensor train. The count may exceed the range that a MATLAB double can represent exactly for large spin systems.

## Physical / mathematical content

The implementation first checks that the input is a `ttclass` object, then multiplies the dimensions returned by `sizes(tt)` using native `int64` arithmetic. It raises an error if the result exceeds `flintmax`, otherwise returns the count as a double.

## Numerical / algorithmic content

The representability check uses MATLAB's `flintmax`; the dimension product is formed with `prod(...,'native')` on `int64` values.

## Parameters / inputs

- tt - tensor train object

## Outputs

- n - an integer-valued double; an error is raised if the count exceeds MATLAB's `flintmax`

## Implementation structure

- Validate that `tt` is a `ttclass` object.
- Compute the product of its sizes in `int64` arithmetic.
- Check against `flintmax` and convert the result to double.
