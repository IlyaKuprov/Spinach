# kernel/conventions/transforms/qter2dcm.m

- Signature: `dcm=qter2dcm(q)`

## Purpose

Converts a quaternion to the 3x3 direction cosine matrix (DCM) used by `euler2dcm.m`.

## Physical / mathematical content

The function first normalizes the quaternion components `(u,i,j,k)`. It then returns `[1-2*(j^2+k^2), 2*(i*j-u*k), 2*(i*k+u*j); 2*(i*j+u*k), 1-2*(i^2+k^2), 2*(j*k-u*i); 2*(i*k-u*j), 2*(j*k+u*i), 1-2*(i^2+j^2)]`. This is the active convention used by `euler2dcm.m`; MATLAB Aerospace Toolbox `quat2dcm()` returns the transpose for the same quaternion.

## Numerical / algorithmic content

The quaternion norm must be at least `sqrt(eps())`; smaller norms cause an error. The DCM acts on a column vector as `v=dcm*v` and transforms a 3x3 tensor as `A=dcm*A*dcm'`.

## Syntax

```matlab
dcm=qter2dcm(q)
```

## Parameters / inputs

- `q` — structure with real numeric scalar fields `u`, `i`, `j`, and `k`.

## Outputs

- `dcm` — 3x3 direction cosine matrix.

## Implementation structure

The function checks the quaternion fields and their values, rejects a norm below `sqrt(eps())`, normalizes the quaternion, and constructs the DCM from its components.
