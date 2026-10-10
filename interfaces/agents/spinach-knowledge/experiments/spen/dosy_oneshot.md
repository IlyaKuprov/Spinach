# experiments/spen/dosy_oneshot.m

- MATLAB source: [experiments/spen/dosy_oneshot.m](https://github.com/IlyaKuprov/Spinach/blob/main/experiments/spen/dosy_oneshot.m)
- Signature: `fid=dosy_oneshot(spin_system,parameters,H,R,K,G,F)`

## Sequence model

The source describes a one-shot DOSY pulse sequence. It forms `L=H+F+1i*R+1i*K`, builds a Y pulse operator for the single specified spin over the spatial Fokker–Planck grid, and starts from the supplied `parameters.rho0`. The implemented pathway is explicit: a 90-degree Y rotation, selection of coherence order -1, a first gradient interval, a 180-degree Y rotation, selection of +1, a second gradient interval, two successive 90-degree Y rotations, and selection of coherence order 0. It then applies the central diffusion interval and refocusing gradient intervals, followed by a 90-degree Y rotation, two more gradient intervals separated by a 180-degree Y rotation, and F2 signal acquisition. No particular prepared state is assumed beyond the caller-supplied `rho0`.

Gradient intervals use `G{1}` and last `g_dur/2` each. The first pair uses the signed terms `+(1+kappa)*g_amp*G{1}` and `-(1-kappa)*g_amp*G{1}`; the central intervals use `-(2*kappa)*g_amp*G{1}`. Each gradient interval is followed by a `g_stab_del` evolution under `L`. The central diffusion evolution lasts `del-4*(g_dur/2)-4*g_stab_del`, with the same central gradient term on its two sides. These are the source's specified intervals; this function does not define a chirped RF pulse.

## Inputs, units, and grids

- `rho0`: caller-prepared initial state; `coil`: detection state; `spins`: a one-element cell array naming the working spin.
- `g_amp`: gradient amplitude, T/m; `g_dur`: gradient pulse width, s; `kappa`: dimensionless bipolar-gradient imbalance, with amplitude ratio `(1+kappa):(1-kappa)`.
- `g_stab_del` and `del`: stabilisation delay and diffusion interval, respectively, in seconds. The code rejects `del < 2*g_dur+4*g_stab_del`.
- `dims`: sample size in m; `npts`: number of spatial discretisation points. `diff` supplies spatially uniform diffusion coefficient/tensor data in m^2/s, or `dxx` supplies voxel-wise diffusion along the sample axis in m^2/s; provide exactly one of these alternatives.
- `npoints`: number of acquired signal points; `sweep`: acquisition sweep width, Hz.
- `H`, `R`, and `K`: Fokker–Planck Hamiltonian, relaxation, and kinetics superoperators; `G`: gradient superoperators; `F`: diffusion and flow superoperator. The function requires the `sphten-liouv` formalism.

## Output axis

The final call uses timestep `1/sweep` and `npoints-1` evolution steps in observable mode. For a single supplied initial state, `fid` is a vector of `npoints` time-domain samples; the spatial Fokker–Planck grid is internal to the calculation, not an extra returned signal axis.

## References

- [Canonical MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/experiments/spen/dosy_oneshot.m)
- [Spinach Wiki: dosy_oneshot.m](https://spindynamics.org/wiki/index.php?title=dosy_oneshot.m)
