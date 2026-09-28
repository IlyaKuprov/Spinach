# kernel/overloads/@polyadic/suffix.m

- Signature: `p=suffix(p,a)`

## Purpose

Associates a suffix factor with a polyadic object. A scalar is multiplied into the final core of each term; a matrix factor is appended to the suffix chain after its dimensions are checked.

## Physical / mathematical content

A polyadic represents a sum of products of core matrices, with optional prefix and suffix factors. The suffix chain records right-side factors of that represented product; it may itself contain a polyadic.

## Numerical / algorithmic content

Scalar factors are absorbed into every term's last core. A non-scalar factor is stored without expanding the polyadic, subject to the inner-dimension check.

## Parameters / inputs

- p -polyadic object
- a -suffix matrix

## Outputs

- p -polyadic object with the factor incorporated as a core scaling or suffix factor
- Note: a suffix can be a polyadic itself.

## Implementation structure

- Checks that p is a polyadic.
- For scalar a, left-multiplies the last core of each term by a.
- Otherwise checks size(p,2) against size(a,1), then appends a to p.suffix.
