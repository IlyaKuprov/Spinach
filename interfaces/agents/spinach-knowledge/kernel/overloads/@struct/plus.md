# kernel/overloads/@struct/plus.m

Direct source: [kernel/overloads/@struct/plus.m](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/overloads/@struct/plus.m)
Wiki: [struct/plus.m](https://spindynamics.org/wiki/index.php?title=struct/plus.m)

- Signature: `str3=plus(str1,str2)`

## Purpose

Add corresponding fields of two structures, recursively applying `+` to nested structures and dispatching non-structure field values to their own overloaded `plus` methods.

## Inputs and output

Both top-level operands must be structures. At each recursive structure level, the implementation checks that both operands have the same number of field names and the same set of names; field ordering need not match. It then evaluates `str1.field + str2.field` for each name. A mismatch, or a non-structure operand at a recursive call, raises `structure topology mismatch.`

`str3` contains the corresponding field results. There is no separate pre-check for leaf dimensions or types: those are governed by the `+` operation selected for each pair of field values.
