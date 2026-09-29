# experiments/imaging/cpmg_dec.m

- Signature: `mri=cpmg_dec(spin_system,parameters,H,R,K,G,F)`
- Canonical MATLAB source: [`experiments/imaging/cpmg_dec.m`](https://github.com/IlyaKuprov/Spinach/blob/main/experiments/imaging/cpmg_dec.m)

## Contract

Call from the `imaging()` context, which supplies `H`, `R`, `K`, `G`, and `F`. The function returns an MRI phantom/image on the spatial grid described by `parameters.npts`; it does not return an acquired FID. The source describes the output as the detection-state amplitude at each sample point.

Required controls are `parameters.dec_time` (positive total free-evolution duration, seconds), `parameters.npulses` (positive integer number of refocusing pi pulses, excluding the initial pi/2), and `parameters.spins` (nonempty cell array of spin labels; only its first entry defines the pulse operator). The input state and grid are supplied through the imaging parameters. The detected state is `parameters.coil_st{1}`; the coil phantom itself is ignored. The source restricts use to the `sphten-liouv` and `zeeman-liouv` formalisms.

## Sequence and propagation

The background generator is `B=H+F+1i*R+1i*K`. Spatially expanded `Lx` and `Ly` pulse operators are formed with the identity across the grid. The sequence applies an ideal pi/2 rotation about `Ly`, then free evolution under `B` for half an interval. With `tau=dec_time/npulses`, it applies `npulses-1` cycles of a pi rotation about `Lx` and a full `tau` evolution, followed by one last pi rotation and a final half-interval. Thus free evolution totals `dec_time`; the documented pulse count excludes the initial pi/2.

The final state is projected with `fpl2phan(rho,parameters.coil_st{1},parameters.npts)`. The returned image follows that spatial grid; no simulated or measured voxel values are asserted here. `G` is checked as a cell array with at least one gradient operator, but is not otherwise used by this function's sequence.

## References

- Spinach Wiki: <https://spindynamics.org/wiki/index.php?title=cpmg_dec.m>
- Source: [`experiments/imaging/cpmg_dec.m`](https://github.com/IlyaKuprov/Spinach/blob/main/experiments/imaging/cpmg_dec.m)
