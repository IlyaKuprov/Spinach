# examples/fundamentals/tensor_structures/polyadic_test_2.m

- Signature: `polyadic_test_2()`

## Purpose

Unit tests for advanced polyadic functionality.

## Physical / mathematical content

- Fundamentals examples. These are unit tests, convention checks, and pedagogical demonstrations of operator algebra, perturbation theory, tensor conventions, symmetry, quadrature, and numerical differentiation.

## Numerical / algorithmic content

- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.

## Implementation structure

- Unit tests for advanced polyadic functionality.
- Get random test matrices
- Build a reference polyadic and its matrix form
- Check constructor, full, inflate, and validate
- Check prefixes, suffixes, size, and emptiness
- Check addition and subtraction paths
- Check multiplication paths
- Check Kronecker products
- Check transpose operations
- Check finiteness and internal non-zero counts
- Check zero-dimension behaviour
- Check nested simplification paths
