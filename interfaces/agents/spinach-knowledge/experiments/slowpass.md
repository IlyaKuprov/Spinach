# experiments/slowpass.m

- MATLAB source: [experiments/slowpass.m](https://github.com/IlyaKuprov/Spinach/blob/main/experiments/slowpass.m)
- Signature: `spectrum=slowpass(spin_system,parameters,H,R,K)`

## What the routine computes

`slowpass` evaluates the response at specified frequency points directly from the Liouvillian resolvent; it does not first calculate a complete FID and Fourier-transform it. The caller supplies the starting state in `parameters.rho0`, the detection state in `parameters.coil`, and the dynamics through `H`, `R`, and `K`.

The routine moves the inputs to the adjoint representation when needed, forms `L=H+1i*R+1i*K`, and obtains subspaces selected for the coil. In each subspace it projects the initial state, coil, and Liouvillian, then evaluates and sums the coil response at each frequency `omega` by solving the shifted Liouvillian system with right-hand side `rho0_subs` and pairing the result with `coil_subs`.

## Identity directions and reaction coupling

In Liouville space, concentration-independent unit vectors identify each substance's identity sector. Spatial contexts supply `parameters.spc_dim`; identities are then embedded in every coordinate of the space-times-spin basis. After subspace projection and normalisation, both directions of coupling between identity and spin order are checked against `spin_system.tols.liouv_zero`.

Only a decoupled identity sector is removed from the initial and detection states and shifted by `-1i*U*U'` in the Liouvillian (one inverse second). This lifts stationary identity poles without changing the spin-order resolvent. If either coupling is nonzero beyond the tolerance, as can occur for selective reactions, the original states and Liouvillian are retained. This is not a general treatment of singular spin-order modes: the relaxation requirement still applies, and coupled stationary modes are not regularised.

Wavefunction inputs do not have a Liouville identity sector and bypass unit-state construction. Hilbert-space density-matrix inputs are converted by `sim2liouv` before this treatment.

## Frequency and output axes

- `parameters.sweep` is a two-element frequency interval in Hz.
- `parameters.npoints` sets the number of points. The code forms `2*pi*linspace(sweep(1),sweep(2),npoints)'`, so the interval is converted to angular frequency for the resolvent.
- The returned `spectrum` is a complex column vector with `npoints` entries, one per equally spaced requested frequency. No time axis or FID is returned.
- At least two frequency points are required. Multiplication by `abs(diff(sweep))*npoints/(npoints-1)` matches the unnormalised FFT amplitude convention.

## Inputs and relaxation condition

- Required state inputs are `parameters.rho0` and `parameters.coil`; `H`, `R`, and `K` are matrices supplied by the context function.
- Relaxation must be present in the dynamics for the source's linear-system calculation, and `R` must not be thermalised. This is a stated input condition, not a report of a run or a measured convergence result.

## References

- [Canonical MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/experiments/slowpass.m)
- [Spinach Wiki: slowpass.m](https://spindynamics.org/wiki/index.php?title=slowpass.m)
