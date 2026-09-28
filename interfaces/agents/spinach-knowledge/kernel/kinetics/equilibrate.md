# kernel/kinetics/equilibrate.m

- Signature: `c=equilibrate(K,c0)`

## Purpose

Computes the equilibrium concentration vector for linear kinetics described by `dc/dt=K*c`.

## Physical / mathematical content

The model conserves total concentration: the columns of `K` must sum to zero, and the returned state has the same total concentration as `c0`.

## Numerical / algorithmic content

The equilibrium satisfies `K*c=0` and `sum(c)=sum(c0)`. Independent reaction components are solved recursively. For each connected component, the rate matrix is rescaled by its largest absolute entry, and the steady-state and mass-conservation equations are assembled into one linear system. The routine checks the system's condition number before solving.

## Parameters / inputs

- `K` — a real square reaction-rate matrix for `dc/dt=K*c`; its column sums must be zero to within `10*eps('double')`.
- `c0` — a non-negative real column vector of initial concentrations, with length equal to the dimension of `K`.

## Outputs

- `c` — the equilibrium concentration vector.

## Implementation structure

After input validation, zero initial concentration is returned directly. Otherwise the code separates independent reactions, normalizes each rate scale, forms the stacked constraints `[ones; K]*c=[sum(c0); zeros]`, and solves if the condition number does not exceed `1/sqrt(eps('double'))`.
