# kernel/overloads/@polyadic/apply.m

[Source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/overloads/@polyadic/apply.m)

`answer=apply(p,x)` evaluates the operator action on the columns of a numeric matrix with `size(p,2)` rows. It returns a numeric block with `size(p,1)` rows and the same column count, even when both the operator and the input are scalar. This distinguishes numerical action from scalar operator scaling.

The method contracts the unopened Kronecker factors and applies stored prefixes and suffixes in matrix-product order. Ordinary multiplication delegates numeric block actions here, including the scalar probes used by `cheap_norm`; no change to the norm-estimation or propagation algorithm is required.
