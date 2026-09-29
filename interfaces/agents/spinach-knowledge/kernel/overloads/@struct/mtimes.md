# kernel/overloads/@struct/mtimes.m

Direct source: [kernel/overloads/@struct/mtimes.m](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/overloads/@struct/mtimes.m)
Wiki: [struct/mtimes.m](https://spindynamics.org/wiki/index.php?title=struct/mtimes.m)

- Signature: `str_out=mtimes(M,str_in)`

## Purpose

Apply left matrix multiplication by `M` to every field of a structure; nested structures are handled recursively through MATLAB's overloaded `mtimes` dispatch.

## Inputs and output

- `M` is numeric. The implementation checks `isnumeric(M)`, not a particular shape.
- `str_in` must be a structure. Its fields are operated on one at a time.
- `str_out` has the corresponding field names, with each field replaced by `M*str_in.field`.

For numeric leaves, ordinary MATLAB matrix-multiplication rules apply: inner dimensions must agree, and the result shape follows those operands. The implementation checks neither each leaf's type nor its dimensions in advance; incompatible fields fail when their `*` operation is evaluated. A nested structure invokes this overload again.
