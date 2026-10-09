# kernel/utilities/sinkhole.m

**Source:** [https://github.com/IlyaKuprov/Spinach/blob/main/kernel/utilities/sinkhole.m](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/utilities/sinkhole.m)

## Purpose

Turns specified states of a spin system into sinkholes: any population reaching them is summed up and stored forever in a frozen state. This is useful for state space restriction diagnostics.

## Behaviour

- Syntax: `L=sinkhole(spin_system,L,states)`.
- The function first runs a consistency check (`grumble`) on the inputs.
- The columns of the Liouvillian corresponding to the sinkhole states are zeroed: `L(:,states)=0`.
- This functionality is only available in the `sphten-liouv` formalism; an error is raised otherwise.
- Consistency enforcement errors if:
  - `spin_system.bas.formalism` is not `'sphten-liouv'`.
  - `L` is not numeric or not a square matrix.
  - The dimension of `L` does not match the dimension of the basis set (`spin_system.bas.offsets(end)`).
  - `states` is not a numeric, real vector of positive integers.
  - Any element of `states` exceeds the state space dimension.

## Inputs and outputs

**Inputs**

- `spin_system` — spin system object.
- `L` — Liouvillian matrix.
- `states` — vector of positive integers specifying the numbers of the states to be set up as sinkholes.

**Outputs**

- `L` — updated Liouvillian matrix with the specified columns zeroed.

## References

- Spinach Wiki: [sinkhole.m](https://spindynamics.org/wiki/index.php?title=sinkhole.m)
