# examples/relaxation_theory/from_md/rwalk_plot.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/relaxation_theory/from_md/rwalk_plot.m`
- Signature: `rwalk_plot()`
- Total lines: 33

## Purpose

A plot of a typical random walk on a sphere.

## Physical / mathematical content

- Relaxation-theory examples. The mathematical backbone is Bloch-Redfield-Wangsness or stochastic Liouville theory, spectral densities, cross-correlation terms, motional models, and extraction of longitudinal/transverse decay behaviour from superoperators.

## Numerical / algorithmic content


## Code-derived implementation details

### Comment-guided execution stages

- Lines 8-9: tau_c and number of steps per tau_c; implemented by `tau_c=6e-9; tc_steps=500`.
- Lines 12-13: Total number of steps; implemented by `tot_nsteps=5000`.
- Lines 15-16: Euler angles of random walk; implemented by `rng('shuffle')`.
- Lines 19-20: Trajectory preallocation; implemented by `traj=zeros(3,size(eulers,1))`.
- Lines 22-23: Trajectory; implemented by `for n=1:size(eulers,1)`.
- Lines 27-28: Plotting; implemented by `plot3(traj(1,:),traj(2,:),traj(3,:))`.

### Control flow inferred from the code

- Line 23: `for` loop over `n=1:size(eulers,1)`.

### Key state/data transformations

- Lines 9: computes `tau_c` using `tau_c=6e-9; tc_steps=500`.
- Lines 10: computes `dt` using `dt=tau_c/tc_steps`.
- Lines 13: computes `tot_nsteps` using `tot_nsteps=5000`.
- Lines 17: computes `eulers` using `eulers=rwalk(tot_nsteps,tau_c,dt)`.
- Lines 20: computes `traj` using `traj=zeros(3,size(eulers,1))`.
- Lines 24: computes `traj(:,n)` using `traj(:,n)=euler2dcm(eulers(n,:))*[0; 0; 1]`.

## Implementation structure

- A plot of a typical random walk on a sphere.
- tau_c and number of steps per tau_c
- Total number of steps
- Euler angles of random walk
- Trajectory preallocation
- Trajectory
- Plotting

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `rng()`, `rwalk()`, `traj()`, `euler2dcm()`, `eulers()`, `plot3()`, `xlim()`, `ylim()`, `zlim()`.
