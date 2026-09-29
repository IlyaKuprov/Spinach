# kernel/utilities/krondelta.m

**Source:** [kernel/utilities/krondelta.m](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/utilities/krondelta.m)

## Purpose

Computes the Kronecker symbol (Kronecker delta) for two integers, returning a logical value indicating whether the two inputs are equal.

## Behaviour

- Syntax: `d=krondelta(a,b)`.
- The function first calls an internal consistency check (`grumble`) on both inputs.
- If the inputs pass validation, the function returns `true` when `a==b` and `false` otherwise.
- The output is created via `true()` or `false()`, so `d` is a logical scalar.

### Input validation

The internal `grumble` subroutine enforces all of the following conditions on both `a` and `b`, and throws an error with the message `'a and b must be real integer scalars.'` if any check fails:

- `isnumeric` must be true.
- `isscalar` must be true.
- `isreal` must be true.
- `mod(a,1)==0` and `mod(b,1)==0` (i.e., the values must be integers).

## Inputs and outputs

**Inputs**

- `a` — an integer number.
- `b` — an integer number.

**Outputs**

- `d` — a logical number (`true` if `a` equals `b`, `false` otherwise).

## References

- Spin Dynamics Wiki page for this function: [krondelta.m](https://spindynamics.org/wiki/index.php?title=krondelta.m)
- The source file closes with the quotation attributed to Leopold Kronecker: "Die ganzen Zahlen hat der liebe Gott gemacht, alles andere ist Menschenwerk."
