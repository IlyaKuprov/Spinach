# kernel/utilities/core_size.m

[Source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/utilities/core_size.m)

`n=core_size(x,dim)` returns the row or column count of a numeric factor or an implicit core construction description. For implicit descriptions it reads `x.dims`; otherwise it uses `size(x,dim)`. The dimension is 1 or 2, and no action is executed.

The helper rejects unsupported core types and requires implicit dimensions to be a finite row of two positive integers. Both the core and requested dimension are validated before accessing metadata.
