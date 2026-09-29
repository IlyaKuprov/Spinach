# kernel/utilities/istraceless.m

## Purpose

Checks whether a matrix is traceless within floating-point precision, returning a logical true or false. Source: [kernel/utilities/istraceless.m](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/utilities/istraceless.m).

## Behaviour

The function calls `grumble(M)` to enforce that the input is numeric, raising the error `'M must be numeric.'` if `~isnumeric(M)` is true. It then computes the working precision as `eps(class(M))`, obtains the cheapest norm of `M` via `cheap_norm(M)`, and returns `A=(abs(trace(M))<=precision*norm_m)`. Thus `A` is true when the absolute value of the trace does not exceed the precision-scaled norm of the matrix.

## Inputs and outputs

**Inputs:**

- `M` — a matrix of any dimension.

**Outputs:**

- `A` — true if the matrix is traceless to appropriate precision, false otherwise.

## References

- Spinach Dynamics Wiki: [istraceless.m](https://spindynamics.org/wiki/index.php?title=istraceless.m)
