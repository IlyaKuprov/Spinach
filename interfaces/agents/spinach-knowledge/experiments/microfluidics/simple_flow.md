# experiments/microfluidics/simple_flow.m

- Signature: `traj=simple_flow(spin_system,parameters,H,R,K,~,F)`

## Purpose

Simple forward evolution experiment for the microfluidics module; trajectory is returned. Syntax: traj=simple_flow(spin_system,parameters,H,R,K,~,F) This sequence must be called from the meshflow() context, which would provide H, R, K, G, and F. Because gradients are not being used, the G input is ignored.

## Physical / mathematical content

- The file uses a Fokker-Planck-style enlarged state space in which spatial or orientational coordinates are promoted to extra dimensions and coupled to spin dynamics through differential operators.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.

## Parameters / inputs

- parameters.npoints -number of points in
- the trajectory
- parameters.rho0 -initial state in Fokker-
- Planck space
- parameters.dt -trajectory time step

## Outputs

- traj -trajectory in the Fokker-Planck space
- Notes: to convert Fokker-Planck space trajectory into R3
- or Liouville space, use fpl2phan and fpl2rho func-
- tions.

## Implementation structure

- Simple forward evolution experiment for the microfluidics
- module; trajectory is returned. Syntax:
- traj=simple_flow(spin_system,parameters,H,R,K,~,F)
- This sequence must be called from the meshflow() context,
- which would provide H, R, K, G, and F. Because gradients
- are not being used, the G input is ignored.
- parameters.npoints -number of points in
- the trajectory
- parameters.rho0 -initial state in Fokker-
- Planck space
- parameters.dt -trajectory time step
- traj -trajectory in the Fokker-Planck space
