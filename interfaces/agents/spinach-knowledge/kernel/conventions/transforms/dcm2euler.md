# kernel/conventions/transforms/dcm2euler.m

- Signatures: `[alpha,beta,gamma]=dcm2euler(dcm)` or `angles=dcm2euler(dcm)`

## Purpose

Extracts ZYZ active Euler angles from a directional cosine matrix. The convention rotates the object rather than the axes; angles are in radians.

## Conversion

Writing `D` for the input matrix, with `Dij` denoting its `(i,j)` entry, the code forms the symmetric Davenport matrix

```text
K = [ D11+D22+D33,  D32-D23,          D13-D31,          D21-D12;
      D32-D23,       D11-D22-D33,     D12+D21,          D13+D31;
      D13-D31,       D12+D21,         D22-D11-D33,     D23+D32;
      D21-D12,       D13+D31,         D23+D32,         D33-D11-D22 ]
```

The eigenvector associated with the largest eigenvalue is taken as the quaternion `q=(q.u,q.i,q.j,q.k)`. The function passes it to `qter2euler`, then sets `alpha=mod(alpha,2*pi)` and `gamma=mod(gamma,2*pi)` (the source comment describes these as wrapped into `[0,2*pi]`). The source describes the resulting angles as those of the proper rotation nearest to the input in the Frobenius norm. It checks the reconstruction with `euler2dcm(alpha,beta,gamma)`; if the input-minus-reconstruction 1-norm exceeds `1e-2`, it displays both matrices and raises an error. Recovering Euler angles from a DCM is generally non-unique: at singular ZYZ rotations, different `alpha`/`gamma` pairs represent the same rotation. The reconstruction check validates the rotation matrix, not recovery of an original angle triple. The source cites I. Y. Bar-Itzhack, J. Guidance Control Dyn. 23 (2000) 1085 for the Davenport-matrix method.

## Inputs and outputs

- `dcm`: real numeric 3-by-3 matrix. Orthogonality and unit determinant deviations above `1e-6` produce warnings; deviations above `1e-2` produce errors. No separate finite-value check is present.
- One requested output: returns `[alpha beta gamma]`, a 1-by-3 row vector in radians. A call with no requested output is also permitted but returns no accessible value.
- Three outputs: the three angles separately in radians. Any other requested output count raises an error.

## References

- [MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/conventions/transforms/dcm2euler.m)
- [Spinach Wiki: dcm2euler.m](https://spindynamics.org/wiki/index.php?title=dcm2euler.m)
