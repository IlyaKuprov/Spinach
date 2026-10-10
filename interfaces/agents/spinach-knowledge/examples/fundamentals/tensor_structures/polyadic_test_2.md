# examples/fundamentals/tensor_structures/polyadic_test_2.m

- Signature: `polyadic_test_2()`
- Source: [`examples/fundamentals/tensor_structures/polyadic_test_2.m`](https://github.com/IlyaKuprov/Spinach/blob/main/examples/fundamentals/tensor_structures/polyadic_test_2.m)

## Purpose

A MATLAB unit test for the `polyadic` matrix representation. It compares the represented operations with explicit dense Kronecker-product matrices; it is not a spin-dynamics simulation.

## Model and representation

The test draws complex random matrices `a` (2×2), `b` (3×3), `c` (2×2 sparse, with `sprandn` density argument 0.75 for each real and imaginary draw), and `d` (3×3). It builds `p=polyadic({{a,b},{c,d}})`, whose reference is `kron(a,b)+kron(c,d)`, a 6×6 matrix. A second object `q` represents `kron(d.',a.')`. No random seed is set.

There is no spin-system specification, basis, Hamiltonian, or physical unit in this example: the factors are generic numerical matrices. The tested tensor construction is a sum of Kronecker products, not a physical-spin tensor basis.

## Use and checks

Run `polyadic_test_2()` with Spinach on the MATLAB path. The no-argument function checks the constructor, `full`, `inflate`, and `validate`; prefix/suffix composition, size and emptiness; addition, subtraction, scalar scaling, and matrix products; Kronecker products; transpose and conjugate transpose; finiteness and nonzero counts; a zero-row matrix; and nested-object simplification. It includes products with dense and sparse operands. Most numerical comparisons use the one-norm error threshold `1e-12`; the GPU comparison uses `1e-10`.

The GPU conversion comparison is conditional on `gpuDeviceCount>0`; otherwise the source prints that this branch was skipped. The source's assertions and messages describe intended checks, not an observed run or pass result.

## Output and limits

The function returns no output argument. It displays section messages when execution reaches them; a failed explicit assertion, validation, or comparison can stop execution with an error. Results depend on unseeded random draws. This compact unit test checks selected operations on small matrices; it does not benchmark performance or establish behaviour for every input, and GPU behaviour is only exercised when a device is available.
