# experiments/microfluidics/simple_flow.m

- MATLAB source: [experiments/microfluidics/simple_flow.m](https://github.com/IlyaKuprov/Spinach/blob/main/experiments/microfluidics/simple_flow.m)
- Spinach Wiki: [simple_flow.m](https://spindynamics.org/wiki/index.php?title=simple_flow.m)

## Purpose

This microfluidics entry point advances an initial state through the supplied Fokker-Planck-space evolution and returns its trajectory. It is called in the meshflow() context; the source describes a parameterised calculation, not a measured flow result.

## Inputs and parameters

Signature: traj=simple_flow(spin_system,parameters,H,R,K,~,F)

- parameters.rho0: required initial state in Fokker-Planck space.
- parameters.dt: required positive finite time step (seconds in Spinach's time convention; the source comment itself only says “time step”).
- parameters.npoints: required positive integer number of trajectory points.
- H, R, K, and F are supplied by the calling context. The routine forms L=H+1i*F+1i*R+1i*K. The placeholder ~ is the gradient argument G, which this source explicitly says is ignored.

## Evolution and output

The routine calls evolution with parameters.rho0, time step parameters.dt, and parameters.npoints-1 evolution steps in trajectory mode. The returned traj is the Fokker-Planck-space state trajectory; it is not an acquired FID. The source points to fpl2phan and fpl2rho for converting such a trajectory to R3 or Liouville space.

This describes the implemented call path; no MATLAB run or measured flow outcome is asserted.
