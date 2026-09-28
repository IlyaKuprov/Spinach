# tests/kernel/test_tensor_vector_suite.m

- Signature: `result=test_tensor_vector_suite()`

## Purpose

Tests tensor, vector, distribution, and relaxation utilities against small analytical cases.

## Physical / mathematical content

- Checks the zero-skew limit of the skew-normal density, isotropic rotational correlation weights and rates, tensor isotropic-shift replacement with anisotropy preserved, and Fokker-Planck vector embedding, projection, and spatial averaging.

## Numerical / algorithmic content

- Checks cubic Hermite interpolation of `x^2` and a skew-normal density with zero skew against the corresponding normal density.
- Checks phantom-to-Fokker-Planck embedding against a Kronecker product, image extraction by observable projection, and spatial averaging back to spin space.
- Checks isotropic rotational correlation weights and rates, a chemical-species state mask, tensor shifts, bidirectional coupling extraction, g-tensor scaling, and an isotropic Zeeman offset.

## Outputs

- `result` - regression test result with explanatory messages.

## Implementation structure

- Creates a regression test result for tensor, vector, and relaxation utilities.
- Compares each utility's output with a small analytical reference using `test_close` or `test_true`.
- Defines minimal spin systems for the correlation-function and tensor-extraction checks.
