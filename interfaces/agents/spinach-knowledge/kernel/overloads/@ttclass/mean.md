# kernel/overloads/@ttclass/mean.m

- Signature: `answer=mean(ttrain,dim)`

## Purpose

Computes the mean of a tensor-train matrix representation along dimension 1 or 2.

## Numerical / algorithmic content

If all mode sizes are one, the function returns the scalar value immediately. When `dim` is omitted, it selects the first non-singleton matrix dimension. Otherwise, it sums the corresponding mode of each core and divides by that mode's size. Only `dim=1` and `dim=2` are accepted; the result is converted to a scalar if all output modes are singleton.

## Parameters / inputs

- ttrain -a tensor train representation of a matrix
- dim -dimension to operate on (dim=1 or dim=2)

## Outputs

- answer -the mean value computed along the specified dimension
