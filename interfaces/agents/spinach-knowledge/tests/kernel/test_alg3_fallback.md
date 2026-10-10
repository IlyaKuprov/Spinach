# tests/kernel/test_alg3_fallback.m

Source: [tests/kernel/test_alg3_fallback.m](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_alg3_fallback.m)

`test_alg3_fallback()` requires a supported GPU and tests the optional custom CUDA gateway boundary using an isolated copy of the production wrapper. It compares all four real/complex operand combinations against native sparse GPU multiplication with no binary, checks argument rejection, and verifies fallback for an invalid binary, and a mocked MATLAB loader error. Mocked storage-layout, allocation, computation, and validation errors must retain their original identifiers rather than trigger fallback.

The test restores the caller's path, working directory, and shadow-warning setting with cleanup guards. It creates only disposable fixtures and never renames, rebuilds, or modifies shipped MEX binaries. Run it separately as documented in `tests/README.md`; it is not registered in the CPU manifest.
