# examples/fundamentals/state_spaces_1.m

- Signature: `state_spaces_1()`

## Purpose

Reproduces Figure 4 of the state-space-restriction accuracy analysis for a pulse-acquire strychnine experiment (https://doi.org/10.1063/1.3624564). The source estimates hours of runtime, much faster on a Tesla A100 GPU.

## Physical / mathematical content

- The spin system is loaded from `strychnine({'1H'})` at 14.1 T. The basis uses the IK-1 approximation, inter-level 7, proximity level 1, scalar-coupling connectivity, and projection +1; the proximity cutoff is 4.0.
- Redfield relaxation is configured with zero equilibrium, `rlx_keep='kite'`, and a 200 ps correlation time. The initial state is `L+` on `1H`.

## Numerical / algorithmic content

- Disables trajectory-level pruning and enables the greedy algorithm; GPU enablement is shown only as a commented option in the source.
- Forms the Liouvillian as the Hamiltonian plus `1i*relaxation`, propagates a trajectory with step 1 ms for 1000 steps, and analyses correlation order.

## Implementation structure

- Creates and bases the system, applies the NMR assumption, constructs the Liouvillian and trajectory, then displays the correlation-order analysis with a logarithmic vertical axis.
