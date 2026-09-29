# experiments/slowpass.m

- MATLAB source: [experiments/slowpass.m](https://github.com/IlyaKuprov/Spinach/blob/main/experiments/slowpass.m)
- Signature: `spectrum=slowpass(spin_system,parameters,H,R,K)`

## What the routine computes

`slowpass` evaluates the response at specified frequency points directly from the Liouvillian resolvent; it does not first calculate a complete FID and Fourier-transform it. The caller supplies the starting state in `parameters.rho0`, the detection state in `parameters.coil`, and the dynamics through `H`, `R`, and `K`.

The routine moves the inputs to the adjoint representation when needed, forms `L=H+1i*R+1i*K`, and obtains subspaces selected for the coil. In each subspace it projects the initial state, coil, and Liouvillian, then evaluates and sums the coil response at each frequency `omega` by solving the shifted Liouvillian system with right-hand side `rho0_subs` and pairing the result with `coil_subs`.

## Frequency and output axes

- `parameters.sweep` is a two-element frequency interval in Hz.
- `parameters.npoints` sets the number of points. The code forms `2*pi*linspace(sweep(1),sweep(2),npoints)'`, so the interval is converted to angular frequency for the resolvent.
- The returned `spectrum` is a complex column vector with `npoints` entries, one per equally spaced requested frequency. No time axis or FID is returned.

## Inputs and relaxation condition

- Required state inputs are `parameters.rho0` and `parameters.coil`; `H`, `R`, and `K` are matrices supplied by the context function.
- Relaxation must be present in the dynamics for the source's linear-system calculation, and `R` must not be thermalised. This is a stated input condition, not a report of a run or a measured convergence result.

## References

- [Canonical MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/experiments/slowpass.m)
- [Spinach Wiki: slowpass.m](https://spindynamics.org/wiki/index.php?title=slowpass.m)
