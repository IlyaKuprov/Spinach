# kernel/overloads/@ttclass/ranks.m

Signature: `ttranks=ranks(ttrain)`

Source: [kernel/overloads/@ttclass/ranks.m](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/overloads/@ttclass/ranks.m) · Wiki: [ttclass/ranks.m](https://spindynamics.org/wiki/index.php?title=ttclass/ranks.m)

For each buffered train, the method returns an `(ncores+1)-by-ntrains` rank array. Column `n` lists the first core's left bond dimension through the last core's left bond dimension, followed by the last core's fourth dimension (the right boundary rank). Thus the first and final entries are the boundary ranks, expected to be 1 for a valid train; this method reads them and does not validate that condition.

It reads core-array dimensions only: it neither evaluates nor materialises tensor entries and leaves the train unchanged. There is no explicit type or shape guard in this method; it directly accesses `ttrain.cores`. The method performs no scalar operation or conjugation.
