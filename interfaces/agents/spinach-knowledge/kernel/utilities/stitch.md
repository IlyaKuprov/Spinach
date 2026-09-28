# kernel/utilities/stitch.m

- Signature: `fid=stitch(spin_system,L,rho_stack,coil_stack,...`

## Purpose

Stitches forward-propagated state trajectories and backward-propagated detection trajectories at the midpoint of a 3D NMR pulse sequence. This obtains a three-dimensional free induction decay from two 2D simulations (http://dx.doi.org/10.1016/j.jmr.2014.04.002).

The complete call is `fid=stitch(spin_system,L,rho_stack,coil_stack,mec_oper,mec_time,t1,t2,t3,tdir)`. If omitted, `tdir` defaults to `'+-'`.

## Parameters / inputs

- `spin_system` — spin-system structure used to construct propagators, apply the propagator cleanup tolerance, select CPU or GPU execution, and report progress.
- `L` — spin-system Liouvillian; a square numerical matrix.
- `rho_stack` — state-vector stack from the forward part of the simulation, with `size(L,1)` rows and `t1.nsteps` columns.
- `coil_stack` — coil-vector stack from the backward part of the simulation, with `size(L,1)` rows and `t3.nsteps` columns.
- `mec_oper` — cell array of midpoint-event operators, each a matrix the same size as `L`; for example, `{Lx,L,Sy}`.
- `mec_time` — cell array of positive real durations for the corresponding midpoint events. It must have the same number of elements as `mec_oper`.
- `t1.nsteps` — positive integer number of time steps in `t1`.
- `t2.nsteps` — positive integer number of time steps in `t2`.
- `t2.timestep` — finite positive real duration of each time step in `t2`.
- `t3.nsteps` — positive integer number of time steps in `t3`.
- `tdir` — optional two-character time-direction specification for state and coil propagation, respectively: `'++'`, `'+-'`, `'-+'`, or `'--'`. The default is `'+-'`.

## Output

- `fid` — three-dimensional free induction decay, indexed as `fid(t3,t2,t1)` and sized `t3.nsteps` by `t2.nsteps` by `t1.nsteps`.

## Algorithm

The function constructs a half-step propagator from `L` using `t2.timestep/2` and builds the midpoint-event propagator by applying the operators in `mec_oper` in sequence for their corresponding durations. At each `t2` step, it contracts the coil and state stacks through the midpoint propagator: `fid(:,k,:) = coil_stack'*(Pm*rho_stack)`. It then advances each stack by a half-step in the directions specified by `tdir`, using the propagator for `'+'` and its conjugate transpose for `'-'`. Stitching runs on the GPU when `'gpu'` is enabled in `spin_system.sys.enable`; otherwise it runs on the CPU. Each computed slice is gathered into `fid`.

Source reference: <https://spindynamics.org/wiki/index.php?title=stitch.m>