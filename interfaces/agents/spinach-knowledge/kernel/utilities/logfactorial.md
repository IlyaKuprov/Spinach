# kernel/utilities/logfactorial.m

- Signature: `lf=logfactorial(n)`

## Purpose

Logarithm of the factorial function. Avoids complications with factorials of large numbers overflowing 64-bit numbers. Syntax: lf=logfactorial(n)

## Physical / mathematical content

- Computes `log(n!)` for non-negative integers without constructing factorials that can overflow double precision.

## Numerical / algorithmic content

- Evaluates MATLAB `gammaln(n+1)` element-wise.

## Parameters / inputs

- n -non-negative integer number

## Outputs

- lf -logarithm of the factorial of n
- Notes: double precision overflow is a persistent problem with
- Clebsch-Gordan coefficients and other objects that in-
- volve factorials.

## Implementation structure

- Validates that every element of `n` is a non-negative integer, then returns `gammaln(n+1)`.
