# kernel/overloads/@polyadic/isempty.m

MATLAB source: <https://github.com/IlyaKuprov/Spinach/blob/main/kernel/overloads/@polyadic/isempty.m>
Spinach Wiki: <https://spindynamics.org/wiki/index.php?title=polyadic/isempty.m>

## Meaning

This overload asks whether the dimensions reported for a polyadic value contain a zero. It evaluates `size(p)` and returns `true` if any returned dimension equals zero; otherwise it returns `false`. Thus the decision is about the object's reported shape, not whether its stored factors contain nonzero numerical values.

The method does not inflate the representation, inspect core or affix values, or perform matrix multiplication. Its shape semantics are delegated to the applicable `size` method. This file defines no broadcasting rule or additional validation of the polyadic object's internal dimensions.

No orientation-specific data or action is handled here.
