# kernel/overloads/@opium/opium.m

- Signature: `M=opium(dim,coeff)`

## Purpose

Object Pretending It is a Unit Matrix (OPIUM). Syntax: M=opium(dim,coeff)

## Physical / mathematical content

## Numerical / algorithmic content

- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.

## Parameters / inputs

- dim -dimension of the unit matrix
- coeff -coefficient in front of the unit matrix

## Outputs

- M -an OPIUM representing the specified matrix

## Implementation structure

- Object Pretending It is a Unit Matrix (OPIUM). Syntax:
- M=opium(dim,coeff)
- dim -dimension of the unit matrix
- coeff -coefficient in front of the unit matrix
- M -an OPIUM representing the specified matrix
- Default properties
- Method description
- Constructor function
- Check consistency
- Store the parameters
- Number of non-zeroes
- Distinguish zero and scaled unit objects
