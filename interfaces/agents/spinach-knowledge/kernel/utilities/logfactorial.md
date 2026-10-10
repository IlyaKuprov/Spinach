# kernel/utilities/logfactorial.m

## Purpose

Computes the natural logarithm of the factorial of a non-negative integer, avoiding direct factorial overflow in 64-bit double precision. The source notes that double precision overflow is a persistent problem with Clebsch-Gordan coefficients and other objects involving factorials.

## Behaviour

- Syntax: `lf=logfactorial(n)`.
- The function validates its input via an internal consistency check (`grumble`) and errors with the message `elements of n must be non-negative integers.` if any element of `n` is non-numeric, non-real, negative, or non-integer.
- The result is computed using MATLAB's built-in `gammaln` function as `lf=gammaln(n+1)`, exploiting the identity that the gamma function generalises the factorial.

## Inputs and outputs

- `n` — non-negative integer number (must be numeric, real, and have all elements non-negative integers).
- `lf` — logarithm of the factorial of `n`.

## References

- Source: <https://github.com/IlyaKuprov/Spinach/blob/main/kernel/utilities/logfactorial.m>
- Spinach Wiki: <https://spindynamics.org/wiki/index.php?title=logfactorial.m>
