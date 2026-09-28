# kernel/overloads/@rcv/plus.m

- Signature: `C=plus(A,B)`

## Purpose

Adds two same-size RCV matrices, or an RCV matrix and a same-size MATLAB sparse matrix. Scalar addition is rejected.

## Mathematical content

For two RCV inputs, the function represents their sum by concatenating their stored row indices, column indices, and values. A MATLAB sparse operand is converted to RCV form before the same addition path is used.

## Numerical / algorithmic content

The operands must have matching dimensions. If either RCV input is marked as GPU-resident, both are converted to GPU arrays before their stored entries are concatenated. A numeric scalar paired with an RCV matrix is explicitly rejected because adding it would make the matrix non-sparse.

## Parameters / inputs

- A -left operand
- B -right operand

## Outputs

- C -sum A+B, RCV sparse matrix

## Implementation structure

- Check operand types; reject scalar-plus-RCV and unsupported combinations.
- For two RCV inputs, check equal dimensions, align GPU residency if needed, then concatenate their stored entries into the result.
- For a MATLAB sparse operand, check dimensions, convert it to RCV, and recurse through the RCV addition path.
