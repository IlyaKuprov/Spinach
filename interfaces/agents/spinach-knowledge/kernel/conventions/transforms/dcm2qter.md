# kernel/conventions/transforms/dcm2qter.m

- Signature: `q=dcm2qter(dcm)`

## Purpose

Converts a direction cosine matrix in the active convention used by `euler2dcm.m` to a unit quaternion. The quaternion components are `q.u` (scalar), `q.i`, `q.j`, and `q.k`.

## Conversion

For `D=dcm`, with `Dij` denoting its `(i,j)` entry, the function forms the four Shepperd pivots

`p = [1+D11+D22+D33; 1+D11-D22-D33; 1-D11+D22-D33; 1-D11-D22+D33]`

and selects the largest. The corresponding component is set to half the square root of that pivot; the remaining components are obtained from the antisymmetric or symmetric off-diagonal combinations below, divided by four times the selected component:

- If `p1` is selected: `q.u=sqrt(p1)/2`, `q.i=(D32-D23)/(4*q.u)`, `q.j=(D13-D31)/(4*q.u)`, `q.k=(D21-D12)/(4*q.u)`.
- If `p2` is selected: `q.i=sqrt(p2)/2`, `q.u=(D32-D23)/(4*q.i)`, `q.j=(D12+D21)/(4*q.i)`, `q.k=(D13+D31)/(4*q.i)`.
- If `p3` is selected: `q.j=sqrt(p3)/2`, `q.u=(D13-D31)/(4*q.j)`, `q.i=(D12+D21)/(4*q.j)`, `q.k=(D23+D32)/(4*q.j)`.
- If `p4` is selected: `q.k=sqrt(p4)/2`, `q.u=(D21-D12)/(4*q.k)`, `q.i=(D13+D31)/(4*q.k)`, `q.j=(D23+D32)/(4*q.k)`.

If `q.u` is negative, all four components are negated; the result is then normalised by its Euclidean norm. This chooses the non-negative-scalar representative of the quaternion double cover.

## Inputs and outputs

- `dcm`: real numeric 3-by-3 matrix. The Frobenius norm of `dcm'*dcm-eye(3)` must not exceed `1e-6`, and `abs(det(dcm)-1)` must not exceed `1e-6`.
- `q`: structure with scalar fields `u`, `i`, `j`, and `k`, normalised to unit length and with `q.u >= 0`.

## References

- [MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/conventions/transforms/dcm2qter.m)
- [Spinach Wiki: dcm2qter.m](https://spindynamics.org/wiki/index.php?title=dcm2qter.m)
