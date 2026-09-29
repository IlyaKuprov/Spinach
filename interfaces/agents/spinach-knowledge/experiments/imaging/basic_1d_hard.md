# experiments/imaging/basic_1d_hard.m

- Signature: `fid=basic_1d_hard(spin_system,parameters,H,R,K,G,F)`
- Context: call from `imaging()`; the context supplies matrices `H`, `R`, `K`, and `F`, plus gradient operators `G`.
- Canonical implementation: [`experiments/imaging/basic_1d_hard.m`](https://github.com/IlyaKuprov/Spinach/blob/main/experiments/imaging/basic_1d_hard.m)

## Sequence

This is a basic one-dimensional hard-pulse imaging sequence. It forms `L=H+F+1i*R+1i*K`, applies a hard 90-degree pulse about y to the initial state, evolves under `L`, applies a hard 180-degree pulse about x, and evolves under `L-ro_grad_amp*G{1}` for prephasing. It then calls `acquire` under `L+ro_grad_amp*G{1}`. Thus the source explicitly uses a negative gradient sign for prephasing and a positive one for readout; the readout gradient amplitude `ro_grad_amp` is in T/m. The pulse operators act on the first entry of `parameters.spins` and are lifted over `parameters.npts` spatial points.

Both evolution calls pass `0.5/sweep` seconds as the evolution step and request `npoints` steps in final-state mode. For example, a sweep of `S` Hz gives a step argument of `0.5/S` seconds; the source separately supplies `npoints`, so that step alone is not a total echo or acquisition duration. `sweep` is a finite positive detection sweep width in Hz, and `npoints` is a finite positive integer FID sample count. These are sequence parameters, not a claim about a reconstructed image or k-space grid.

## Inputs and output

Use one of the supported formalisms `sphten-liouv` or `zeeman-liouv`. `H`, `R`, `K`, and `F` must be same-size matrices; `G` must be a cell array with at least one gradient operator; and `parameters.rho0` must be numeric. `parameters.spins` is a nonempty cell array of character strings, with the first entry used for the pulse. The routine checks positive-integer `npts` (the spatial lifting count) and `npoints`, as well as positive real `sweep` and real scalar `ro_grad_amp`.

The source header lists `parameters.offset` in Hz. The `imaging()` context applies `frqoffset(spin_system,H,parameters)` before passing `H` to this sequence, so the offset is already in the Hamiltonian. The sequence does not read the offset directly; `F` is the separate diffusion/flow generator from `v2fplanck`, combined with `H` in `L`. The return value is the FID produced by `acquire` for the requested `npoints`; this function does not Fourier transform it or return an image tensor.

## References

- Spin Dynamics Wiki: [`basic_1d_hard.m`](https://spindynamics.org/wiki/index.php?title=basic_1d_hard.m).
- MATLAB implementation: [`experiments/imaging/basic_1d_hard.m`](https://github.com/IlyaKuprov/Spinach/blob/main/experiments/imaging/basic_1d_hard.m).
