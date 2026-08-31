# examples/microfluidics/plain_reaction.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/microfluidics/plain_reaction.m`
- Signature: `plain_reaction()`
- Total lines: 52

## Purpose

Non-linear reaction kinetics in a situation when there is no hydrodynamics, diffusion, or spin dynamics. This is in- tended as a stepping stone to the more complicated cases in the same directory of the Spinach example set. Calculation time: seconds.

## Physical / mathematical content

- Microfluidics examples. The coupled model is spin dynamics plus advection-diffusion-reaction transport on a mesh or regular grid. Numerical issues include finite-difference operators, mesh interpolation, and coupled reaction-flow evolution.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 13-14: No spin system here; implemented by `spin_system=bootstrap()`.
- Lines 16-17: Rate constants, mol/(L*s); implemented by `k1=0.5`.
- Lines 20-21: Cycloaddition reaction generator, including solvent; implemented by `K=@(t,x)(1i*[-k1*x(2)-k2*x(2) 0 0 0 0`.
- Lines 27-28: Kinetic time grid, 20 seconds; implemented by `nsteps=200; tmax=20; dt=tmax/nsteps`.
- Lines 31-32: Preallocate concentration trajectory; implemented by `x=zeros(5,nsteps+1)`.
- Lines 34-35: Initial concentrations, mol/L; implemented by `x(:,1)=[0.6; 0.5; 0.0; 0.0; 18.1]`.
- Lines 37-38: Concentration dynamics; implemented by `for n=1:nsteps`.
- Lines 42-43: Plot concentrations, excluding solvent; implemented by `kfigure(); plot(time_axis,real(x(1:4,:)))`.

### Control flow inferred from the code

- Line 38: `for` loop over `n=1:nsteps`.

### Key state/data transformations

- Lines 14: computes `spin_system` using `spin_system=bootstrap()`.
- Lines 17: computes `k1` using `k1=0.5`.
- Lines 18: computes `k2` using `k2=0.1`.
- Lines 21: computes `K` using `K=@(t,x)(1i*[-k1*x(2)-k2*x(2) 0 0 0 0`.
- Lines 28: computes `nsteps` using `nsteps=200; tmax=20; dt=tmax/nsteps`.
- Lines 29: computes `time_axis` using `time_axis=linspace(0,tmax,nsteps+1)`.
- Lines 32: computes `x` using `x=zeros(5,nsteps+1)`.
- Lines 35: computes `x(:,1)` using `x(:,1)=[0.6; 0.5; 0.0; 0.0; 18.1]`.
- Lines 39: computes `x(:,n+1)` using `x(:,n+1)=step(spin_system,{K,n*dt,'LG4'},x(:,n),dt)`.

## Implementation structure

- Non-linear reaction kinetics in a situation when there is
- no hydrodynamics, diffusion, or spin dynamics. This is in-
- tended as a stepping stone to the more complicated cases
- in the same directory of the Spinach example set.
- Calculation time: seconds.
- No spin system here
- Rate constants, mol/(L*s)
- Cycloaddition reaction generator, including solvent
- Kinetic time grid, 20 seconds
- Preallocate concentration trajectory
- Initial concentrations, mol/L
- Concentration dynamics

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `bootstrap()`, `step()`, `kfigure()`, `kxlabel()`, `kylabel()`, `klegend()`.
