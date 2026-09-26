# kernel/utilities/stitch.m

- Signature: `fid=stitch(spin_system,L,rho_stack,coil_stack,...`

## Purpose

Stitching function for bidirectionally propagated 3D NMR pulse sequences. Propagate your initial condition forward to some mid- point, propagate your detection state backward to the same mid- point and use this function to obtain the 3D free induction de- cay (http://dx.doi.org/10.1016/j.jmr.2014.04.002) at the price of two 2D simulations. Syntax: fid=stitch(spin_system,L,rho_stack,coil_stack,... mtp_oper,mtp_time,t1

## Physical / mathematical content

- General mathematical and infrastructure utilities. This area contains finite differences, perturbation theory, graph algorithms, spectral densities, tensor algebra, hash/report helpers, and other reusable numerical components.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.

## Parameters / inputs

- L -spin system Liouvillian
- rho_stack -state vector stack from the forward part of
- the simulation
- coil_stack -coil vector stack from the backward part of
- the sumulation
- mec_oper -cell array of operators in the midpoint event
- chain, e.g. {Lx,L,Sy}
- mec_time -cell array of durations of each event at the
- midpoint of the t2 evolution period
- t1.nsteps -number of time steps in t1
- t2.nsteps -number of time steps in t2
- t2.timestep -duration of each time step in t2
- t3.nsteps -number of time steps in t3
- tdir -time direction for state and coil propagation,
- the default is '+-'

## Outputs

- fid -three-dimensional free induction decay

## Implementation structure

- Stitching function for bidirectionally propagated 3D NMR pulse
- sequences. Propagate your initial condition forward to some mid-
- point, propagate your detection state backward to the same mid-
- point and use this function to obtain the 3D free induction de-
- cay (http://dx.doi.org/10.1016/j.jmr.2014.04.002) at the price
- of two 2D simulations. Syntax:
- fid=stitch(spin_system,L,rho_stack,coil_stack,...
- mtp_oper,mtp_time,t1,t2,t3)
- L -spin system Liouvillian
- rho_stack -state vector stack from the forward part of
- the simulation
- coil_stack -coil vector stack from the backward part of
