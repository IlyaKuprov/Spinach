# kernel/utilities/stitch.m

**Source:** [kernel/utilities/stitch.m](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/utilities/stitch.m)

## Purpose

Stitches forward- and backward-propagated trajectories to obtain a three-dimensional free induction decay (FID) for 3D NMR pulse sequences at the cost of two 2D simulations. The initial condition is propagated forward to a midpoint, the detection state is propagated backward to the same midpoint, and this function combines the two stacks.

## Behaviour

- Syntax: `fid=stitch(spin_system,L,rho_stack,coil_stack,mec_oper,mec_time,t1,t2,t3,tdir)`.
- Sets the default time direction `tdir` to `'+-'` when the argument is absent.
- Runs a consistency check (`grumble`) on all inputs before computation.
- Preallocates `fid` as a complex array of size `t3.nsteps x t2.nsteps x t1.nsteps`.
- Computes a half-step propagator `P = propagator(spin_system,L,t2.timestep/2)` and its conjugate transpose `Pct = P'`.
- Builds the midpoint event propagator `Pm` as the ordered product of propagators for each event in `mec_oper` with durations from `mec_time`, starting from a sparse identity of the size of `L`, then applies `clean_up` with `spin_system.tols.prop_chop`.
- If `'gpu'` is listed in `spin_system.sys.enable`, uploads `P`, `Pm`, `Pct`, `rho_stack` and `coil_stack` to the GPU via `gpuArray` and reports that stitching will be done on GPU; otherwise reports CPU stitching.
- For each of the `t2.nsteps` steps, reports progress, computes `fid(:,k,:) = gather(coil_stack'*(Pm*rho_stack))`, then advances both stacks according to `tdir`:
  - `'++'`: forward evolution of both `rho_stack` and `coil_stack` using `P`.
  - `'+-'`: forward evolution of `rho_stack` with `P`, backward evolution of `coil_stack` with `Pct`.
  - `'-+''`: backward evolution of `rho_stack` with `Pct`, forward evolution of `coil_stack` with `P`.
  - `'--'`: backward evolution of both stacks using `Pct`.
  - Any other value raises the error `'invalid time direction specification in tdir'`.

## Inputs and outputs

Inputs:

- `spin_system` — Spinach spin system object.
- `L` — spin system Liouvillian; must be a square matrix.
- `rho_stack` — state vector stack from the forward part of the simulation; numerical array of dimensions `size(L,1) x t1.nsteps`.
- `coil_stack` — coil vector stack from the backward part of the simulation; numerical array of dimensions `size(L,1) x t3.nsteps`.
- `mec_oper` — cell array of operators in the midpoint event chain (e.g. `{Lx,L,Sy}`); each element must be a square matrix of the same dimensions as `L`.
- `mec_time` — cell array of durations of each event at the midpoint of the t2 evolution period; each element must be a positive real scalar, and the count must match `mec_oper`.
- `t1` — struct with field `nsteps`, a positive real integer giving the number of time steps in t1.
- `t2` — struct with fields `nsteps` (positive real integer, number of time steps in t2) and `timestep` (finite positive real scalar, duration of each time step in t2).
- `t3` — struct with field `nsteps`, a positive real integer giving the number of time steps in t3.
- `tdir` — optional time direction for state and coil propagation; one of `'++'`, `'+-'`, `'-+'`, `'--'`, default `'+-'`.

Output:

- `fid` — three-dimensional free induction decay, size `t3.nsteps x t2.nsteps x t1.nsteps`.

## References

- Spinach Wiki page: [stitch.m](https://spindynamics.org/wiki/index.php?title=stitch.m)
- Method reference: [http://dx.doi.org/10.1016/j.jmr.2014.04.002](http://dx.doi.org/10.1016/j.jmr.2014.04.002)
