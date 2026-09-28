# kernel/utilities/sinkhole.m

- Signature: `L=sinkhole(spin_system,L,states)`

## Purpose

Turns the specified states into sinkholes: population reaching them is summed up and stored forever in a frozen state. This is useful for state space restriction diagnostics.

## Numerical / algorithmic content

- Checks input consistency, then sets the columns of `L` corresponding to `states` to zero.

## Parameters / inputs

- `spin_system` — spin system; its basis and formalism are used for consistency checks.
- `L` — Liouvillian matrix; must be square and match the dimension of the basis set.
- `states` — vector of positive integers specifying the states to be set up as sinkholes; indices must not exceed the state space dimension.

## Output

- `L` — updated Liouvillian matrix.

## Note

This functionality is only available in `sphten-liouv` formalism.

[Source documentation](https://spindynamics.org/wiki/index.php?title=sinkhole.m)