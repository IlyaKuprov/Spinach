# kernel/utilities/logfactorial.m

- Signature: `lf=logfactorial(n)`

## Purpose

Logarithm of the factorial function. Avoids complications with factorials of large numbers overflowing 64-bit numbers. Syntax: lf=logfactorial(n)

## Physical / mathematical content

- General mathematical and infrastructure utilities. This area contains finite differences, perturbation theory, graph algorithms, spectral densities, tensor algebra, hash/report helpers, and other reusable numerical components.

## Numerical / algorithmic content

## Parameters / inputs

- n -non-negative integer number

## Outputs

- lf -logarithm of the factorial of n
- Notes: double precision overflow is a persistent problem with
- Clebsch-Gordan coefficients and other objects that in-
- volve factorials.

## Implementation structure

- Logarithm of the factorial function. Avoids complications with
- factorials of large numbers overflowing 64-bit numbers. Syntax:
- lf=logfactorial(n)
- n -non-negative integer number
- lf -logarithm of the factorial of n
- Clebsch-Gordan coefficients and other objects that in-
- volve factorials.
- Check consistency
- Use the built-in log(gamma(n)) function
- Consistency enforcement
- For ten years or so, my name was "that jerk". But
- that was a promotion. Before, I was "Who's he?"
