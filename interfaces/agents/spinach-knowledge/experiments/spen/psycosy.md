# experiments/spen/psycosy.m

- Signature: `fid=psycosy(spin_system,parameters,H,R,K,G,F)`
- Source: [Kenwright, spatially encoded COSY, Figure 4](https://doi.org/10.1002/mrc.4727); [Spinach documentation](https://spindynamics.org/wiki/index.php?title=psycosy.m).

## Purpose and output

Implements Alan Kenwright's spatially encoded COSY sequence. Returns `fid`, a two-dimensional free induction decay.

## Inputs

- `spin_system`: Spinach spin system; the implementation requires `sphten-liouv` formalism.
- `parameters.sweep`: sweep width, Hz.
- `parameters.npoints`: two point counts, one for each dimension.
- `parameters.spins`: nuclei on which the sequence runs, specified as `{'1H'}`, `{'13C'}`, etc.; the implementation requires a one-element cell array.
- `parameters.tmix`: mixing time, s.
- `parameters.gamp`: gradient amplitude, T/m.
- `parameters.sal_ang`: saltire chirp flip angle, degrees.
- `parameters.sal_dur`: saltire chirp pulse width, s.
- `parameters.sal_del`: chirp-pulse gradient duration, s.
- `parameters.sal_swp`: saltire chirp sweep width, Hz.
- `parameters.sal_npt`: number of points in the saltire chirp.
- `parameters.sal_smf`: saltire chirp smoothing factor.
- `parameters.rho0`: initial state vector.
- `parameters.coil`: detection state vector.
- `H`: Fokker–Planck Hamiltonian from the imaging context.
- `R`: Fokker–Planck relaxation superoperator from the imaging context.
- `K`: Fokker–Planck kinetics superoperator from the imaging context.
- `G`: Fokker–Planck gradient superoperators from the imaging context; a three-element cell array, with `G{1}` used in this sequence.
- `F`: Fokker–Planck diffusion and flow superoperator from the context.

## Sequence and computation

1. Check input consistency. Set the coherent-evolution timestep `parameters.delta=1/(4*parameters.sweep)` and compose `L=H+F+1i*R+1i*K`.
2. Construct `Lp=operator(spin_system,'L+',parameters.spins{1})`, then spatially extend the pulse operators as `Lx=kron(speye(prod(parameters.npts)),(Lp+Lp')/2)` and `Ly=kron(speye(prod(parameters.npts)),(Lp-Lp')/2i)`.
3. Generate `[Cx,Cy]=chirp_pulse(parameters.sal_npt,parameters.sal_dur,parameters.sal_swp,parameters.sal_smf,'saltire')`. Normalize both components by `norm_factor=max(Cx)`. Calculate `q_beta=-(2*log(cosd(parameters.sal_ang)/2+1/2))/pi` and the saltire RF field strength in Hz, `rfbeta=sqrt(parameters.sal_swp*q_beta/(2*pi*parameters.sal_dur))`; calibrate both components by multiplying by `2*pi*rfbeta`.
4. Apply a hard `pi/2` pulse about `Lx` to `parameters.rho0`. Evolve the first half of `t1` under `L` with timestep `0.5/parameters.sweep`, `parameters.npoints(1)-1` steps and `'trajectory'` mode. Select `+1` coherence on `parameters.spins{1}`, apply a hard `pi` pulse about `Lx`, then select `-1` coherence.
5. Set each chirp segment duration to `parameters.sal_dur/numel(Cx)`. Apply the **first PSYCHE chirp** with operators `{Lx,Ly}` and phases `{Cx,+Cy}` under `L+parameters.gamp*G{1}`, using `'expv-pwc'`. Evolve under that same gradient-containing generator for `parameters.sal_del-2*parameters.sal_dur` in one `'final'` step, then select `0` coherence. Apply the **second PSYCHE chirp** with `{Lx,Ly}` and phases `{Cx,-Cy}`, again under `L+parameters.gamp*G{1}` using the same segment durations and `'expv-pwc'`.
6. Evolve the second half of `t1` under `L` with timestep `0.5/parameters.sweep`, `parameters.npoints(1)-1` steps and `'refocus'` mode. Select `+1` coherence; evolve under the full `L` for `parameters.tmix` in one `'final'` step; select `+1` coherence again. Apply a final hard `pi/2` pulse about `Lx`.
7. Acquire the F2 signal by evolving under `L` with `parameters.coil`, timestep `1/parameters.sweep`, `parameters.npoints(2)-1` steps and `'observable'` mode.

## Consistency requirements

`H`, `R`, `K` and `F` must be numeric, equally sized matrices; `G` must contain three numeric gradient matrices. `parameters.rho0` and `parameters.coil` must be numeric column vectors with the same row count as `H`. `parameters.gamp` and `parameters.sal_ang` must be finite real scalars; `parameters.sweep`, `parameters.sal_dur` and `parameters.sal_swp` must be positive real scalars. `parameters.npoints` must contain two integers greater than 1; `parameters.tmix` must be non-negative; `parameters.sal_del` must exceed `2*parameters.sal_dur`; `parameters.sal_npt` must be a positive integer; and `parameters.sal_smf` must be between 0 and 50.