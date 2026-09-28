# kernel/optimcon/aux_mat.m

- Signature: `[auxm_l,auxm_r]=aux_mat(drifts,controls,cc_comm_idx,...`

## Purpose

Build the left- and right-edge auxiliary matrices used to differentiate the trapezium-product quadrature propagator

`expm(-1i*((HL+HR)/2+(1i*dt/12)*[HL,HR])*dt)`

with respect to control coefficients in the interval generators. The derivative formula follows Eq. 16 of Goodwin and Kuprov ([doi:10.1063/1.4928978](https://doi.org/10.1063/1.4928978)).

## Syntax

```matlab
[auxm_l,auxm_r]=aux_mat(drifts,controls,cc_comm_idx,...
                       cc_comm,dt,cL,cR,k,j)
```

## Inputs

- `drifts` — two drift-generator matrices, ordered left then right.
- `controls` — cell array of `K` control-generator matrices.
- `cc_comm_idx` — `K`-by-`K` logical matrix indicating which control commutators are nonzero.
- `cc_comm` — `K`-by-`K` cell array of control commutation relations.
- `dt` — interval duration in seconds.
- `cL`, `cR` — control coefficients at the left and right interval edges.
- `k` — index of the control generator being differentiated.
- `j` — optional index of a second control generator; supply it for the mixed-derivative 3-by-3 auxiliary matrices.

## Outputs

- `auxm_l`, `auxm_r` — left- and right-edge auxiliary matrices for the requested propagator derivative.

## Construction

The routine forms the interval generator from the average of the two edge generators and their commutator term. With no second index `j`, it returns 2-by-2 block matrices; when `j` is supplied, it returns 3-by-3 block matrices for the mixed derivative. The left and right matrices use the corresponding edge-specific directional derivatives.
