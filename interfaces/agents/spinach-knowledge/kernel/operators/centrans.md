# kernel/operators/centrans.m

Direct source: [kernel/operators/centrans.m](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/operators/centrans.m)

- Signature: `A=centrans(mult,type)`

## Purpose

Construct a sparse complex `mult-by-mult` matrix whose only populated entries are on the two central basis indices. The routine uses `r=mult/2` and `s=r+1` (MATLAB's 1-based indexing) and accepts even integer `mult>=2`. The accepted character values are `x`, `y`, `z`, `+`, and `-`.

## Matrix entries and normalisation

The returned entries are exactly:

| `type` | Nonzero matrix entries |
| --- | --- |
| `x` | `A(r,s)=A(s,r)=0.5` |
| `y` | `A(r,s)=-0.5i` and `A(s,r)=+0.5i` |
| `z` | `A(r,r)=+0.5` and `A(s,s)=-0.5` |
| `+` | `A(r,s)=1` |
| `-` | `A(s,r)=1` |

No additional scale factor is applied; the matrix is converted to complex form before return. The source describes these as central-transition operators for half-integer spins in the Pauli basis.

## Operator action

The function returns an operator matrix only. It does not construct left or right multiplication, a superoperator, or a propagator; any action is determined by how a caller uses the matrix.

## Reference

- [Spin Dynamics documentation for `centrans.m`](https://spindynamics.org/wiki/index.php?title=centrans.m)
