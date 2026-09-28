# experiments/esr_dipolar/deer_analyt.m

- Signature: `deer=deer_analyt(D,J,t)`

## Purpose

Analytical expression for a DEER trace for two spins in the presence of dipolar and exchange coupling.

## Parameters / inputs

- `D` — dipolar coupling in angular frequency units; the coefficient in front of `(1-3*cos(theta)^2)*Lz*Sz` in the spin Hamiltonian. Must be a positive real scalar.
- `J` — exchange coupling in angular frequency units, using the NMR convention (no factor of 2 in front); the coefficient in front of `L*S` in the spin Hamiltonian. Must be a real scalar.
- `t` — array of time points in seconds. Must contain non-negative real numbers.

## Output

- `deer` — array of DEER form factor values with the same dimensions as `t`. At `t==0`, its value is 1.

## Implementation structure

The trace is calculated using [Kuprov's formula](http://dx.doi.org/10.1038/ncomms14842), with `fresnelc` and `fresnels`. The zero-time indeterminacy is removed by setting `deer(t==0)=1`.

<https://spindynamics.org/wiki/index.php?title=deer_analyt.m>