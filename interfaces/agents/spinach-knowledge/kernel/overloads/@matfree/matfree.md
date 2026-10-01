# kernel/overloads/@matfree/matfree.m

[MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/overloads/@matfree/matfree.m)

`core=matfree(dims,forward,adjoint,real_flag)` wraps a matrix-free linear action as a factor in a polyadic Kronecker product. `dims=[rows columns]` contains finite integer dimensions, each at least two. `forward` accepts full `columns`-by-`n` blocks and returns `rows`-by-`n` blocks; `adjoint` has the reversed shapes and must be its Hermitian adjoint. `real_flag` is a logical scalar describing whether the represented matrix is real.

The actions must be linear, finite, and device-preserving. The wrapper validates the handles and metadata, not their captured data or numerical behaviour. `allfinite` checks the stored scaling coefficient; correctness of the actions remains the caller's responsibility. Scalar multiplication scales that coefficient. Transpose and adjoint swap the actions appropriately. The structural `nnz` result is a zero/nonzero marker, not a count of matrix entries. No identity is inferred from an opaque handle.

The existing polyadic contraction applies cores to reshaped blocks. Implicit compositions are buffered instead of eagerly multiplying cores, and singleton implicit cores retain the polyadic norm-estimation route. `full` and `inflate` deliberately reject materialisation. GPU dispatch forwards GPU blocks to the supplied actions without inspecting or transferring their closures. The FFT rotor derivative is an application of this interface.
