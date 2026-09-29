# kernel/optimcon/aux_mat.m

[MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/optimcon/aux_mat.m) · [Spin Dynamics Wiki](https://spindynamics.org/wiki/index.php?title=aux_mat.m)

## Purpose

Build the left- and right-edge block matrices used to calculate directional derivatives of the trapezium-product quadrature propagator. This helper constructs derivative-support matrices; it does not evaluate an optimisation objective or impose constraints. The interval propagator is formed from

The left-edge and right-edge generators as `expm(-1i*((HL+HR)/2+(1i*dt/12)*[HL,HR])*dt)`. This is the formula described by Eq. 16 of Goodwin and Kuprov ([DOI: 10.1063/1.4928978](https://doi.org/10.1063/1.4928978)).

## Syntax

`[auxm_l,auxm_r]=aux_mat(drifts,controls,cc_comm_idx,cc_comm,dt,cL,cR,k,j)`

## Inputs

- `drifts` — cell array containing exactly two square numeric generator matrices, for the left and right interval edges in that order. They must have the same dimensions.
- `controls` — cell array of `K` square numeric control-generator matrices, each with the drift-matrix dimensions.
- `cc_comm_idx` — logical `K`-by-`K` commutation mask used with `cc_comm` to gate control-control terms.
- `cc_comm` — `K`-by-`K` cell array of numeric square commutation matrices, each with the drift-matrix dimensions.
- `dt` — positive real numeric scalar interval duration, in seconds.
- `cL` and `cR` — real numeric arrays with one coefficient per control, at the left and right edges. The source checks their element counts, not a particular vector orientation.
- `k` — positive integer control index from `1` through `K`; differentiation is with respect to this control.
- `j` — optional second control index. It must be an integer from `0` through `K`. If omitted, it defaults to `0`.

## Outputs

- `auxm_l` and `auxm_r` — left- and right-edge auxiliary matrices. For an `M`-by-`M` generator, each is `2M`-by-`2M` when `j=0`, and `3M`-by-`3M` when a nonzero `j` is supplied.

## Construction

The code combines each drift with its edge-specific control coefficients, then forms `G=(G_L+G_R)/2 + 1i*dt*(G_L*G_R-G_R*G_L)/12`. The selected control's left- and right-edge directional generator derivatives occupy the first upper off-diagonal block in the corresponding output. With `j=0`, the diagonal contains three copies of `G` and the lower block is zero. With a second index, a 3-by-3 block arrangement also carries the second control's edge derivative for mixed-derivative calculations. The commutation index and cell array provide the control-control terms used in those derivatives.

Only `dt` has a unit stated by the source (seconds); the source does not prescribe independent units for control coefficients or generator matrices.
