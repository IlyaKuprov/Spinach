# kernel/pulses/iserstep.m

- Signature: `rho_b=iserstep(spin_system,LTM,rho_a,dt)`

## Purpose

Lie-group and Runge-Kutta-Munthe-Kaas solvers for the Lie equa- tion. LG methods are implementations of Equation A.1, with mi- nor typos fixed, from The key difference from step() function is that the Liouvillian can depend on the density matrix. Syntax: rho_b=iserstep(spin_system,{L,t,method},rho_a,dt)

## Physical / mathematical content

- Pulse and waveform utilities. These files encode shaped RF pulses, gradient events, rotating-frame transformations, resonator response, and Lie-group integration of time-dependent driven dynamics.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.

## Parameters / inputs

- spin_system -Spinach data structure from create.m
- and basis.m constructors
- L -a handle to a function L(t,rho) that must take
- time and state vector, and return the evolution
- generator (in rad/s) of the Lie equation:
- d_rho/d_t = -i*L(t,rho)*rho
- rho_a -state vector at the start of the evolution
- period
- t -time at the start of the evolution, seconds
- dt -evolution time step, seconds
- method -'PWCL', 'PWCM', 'RKMK4', 'RKMK-DP5',
- 'RKMK-DP8', or 'LG4'; the latter one
- has a good balance of efficiency and
- numerical accuracy

## Outputs

- rho_b -state vector at the end of the evolution
- time step

## Implementation structure

- Lie-group and Runge-Kutta-Munthe-Kaas solvers for the Lie equa-
- tion. LG methods are implementations of Equation A.1, with mi-
- nor typos fixed, from
- The key difference from step() function is that the Liouvillian
- can depend on the density matrix. Syntax:
- rho_b=iserstep(spin_system,{L,t,method},rho_a,dt)
- spin_system -Spinach data structure from create.m
- and basis.m constructors
- L -a handle to a function L(t,rho) that must take
- time and state vector, and return the evolution
- generator (in rad/s) of the Lie equation:
- d_rho/d_t = -i*L(t,rho)*rho
