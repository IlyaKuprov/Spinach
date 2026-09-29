# kernel/overloads/@opium/size.m

- Signature: `varargout=size(op,dim)`

## Meaning and behaviour

The represented matrix has shape `[op.dim op.dim]`. With no `dim` argument, one output is that two-element shape vector; with two outputs, each output is `op.dim`. With a requested dimension, the method returns `op.dim` for dimension 1 or 2. The implementation checks that a supplied `dim` is scalar and a member of `[1 2]`; other values raise the source's polyadic-dimension error. Other call/output forms reach `invalid call syntax.`.

This query reports the represented square shape without constructing the matrix. The source supports only the stated two dimensions; it does not define higher-dimension, cell-array, or broadcasting behaviour.

## Input and output

- Input: `op`, an `opium` object; optional scalar `dim` equal to 1 or 2.
- Output: the shape vector, both separate dimensions, or the requested single dimension, as described above.

## Source links

- MATLAB source: [kernel/overloads/@opium/size.m](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/overloads/@opium/size.m)
- Existing Wiki page: [opium/size.m](https://spindynamics.org/wiki/index.php?title=opium/size.m)
