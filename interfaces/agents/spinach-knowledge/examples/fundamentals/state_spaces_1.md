# examples/fundamentals/state_spaces_1.m

- Signature: `state_spaces_1()`
- Source: [`examples/fundamentals/state_spaces_1.m`](https://github.com/IlyaKuprov/Spinach/blob/main/examples/fundamentals/state_spaces_1.m)
- Reference: [state-space restriction accuracy analysis, Figure 4](https://doi.org/10.1063/1.3624564)

## Model and question

This pulse-acquire strychnine example tracks the density operator by spin-correlation order in a restricted basis. Its question is how the represented state-space correlation content evolves under the configured Hamiltonian and relaxation, not whether an angular quadrature is converged. The source describes the calculation as a reproduction of Figure 4; that is the example's stated aim, not a result independently established here.

The proton system comes from `strychnine({'1H'})` at 14.1 T. The basis is `sphten-liouv` with `IK-1`, inter-level 7, proximity level 1, scalar-coupling connectivity, and projection `{1}`; the proximity cutoff is 4.0. The algorithm options disable `trajlevel` and enable `greedy`. Redfield relaxation is configured with zero equilibrium, `rlx_keep='kite'`, and a 200 ps correlation time. The initial state is proton `L+`; the Hamiltonian is built under the NMR assumption.

## Propagation and output

The source forms `L=hamiltonian(spin_system)+1i*relaxation(spin_system)` and calls trajectory-mode `evolution` with a 1 ms step for 1000 steps. `trajan(...,'correlation_order')` plots correlation-order content; the displayed axes are 0 to 1000 on x and 1e-7 to 10 on a logarithmic y scale. These are display limits, not declared pass/fail thresholds.

The source estimates hours of runtime and says a Tesla A100 is much faster; GPU enablement is only a commented option. It supplies no numerical acceptance tolerance or quadrature rule. No assertion about agreement with Figure 4 or successful execution follows from the source alone.
