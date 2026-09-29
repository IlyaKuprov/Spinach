# kernel/overloads/@double/inflate.m

- MATLAB implementation: [kernel/overloads/@double/inflate.m](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/overloads/@double/inflate.m)

## Purpose and use

B=inflate(A) is the @double overload of Spinach's inflate interface. It deliberately does not inflate, reshape, or otherwise transform an ordinary double array: the output B is the input A unchanged. This dummy overload lets code using the polyadic inflate command also call it for proper numerical arrays without changing them.

## Inputs and outputs

- A: the double-array argument supplied to this overload; no units or shape restrictions are imposed in this file.
- B: the same value returned unchanged.

## Limitation

This implementation is intentionally a no-op; do not expect it to add dimensions or convert an array into a polyadic representation.

Source documentation: [double/inflate.m](https://spindynamics.org/wiki/index.php?title=double/inflate.m).
