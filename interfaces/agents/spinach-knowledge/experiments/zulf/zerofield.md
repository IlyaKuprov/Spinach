# experiments/zulf/zerofield.m

- Signature: `fid=zerofield(spin_system,parameters,H,R,K)`

## Purpose

Simulates a gamma-weighted pulse-and-acquire experiment in zero field. The initial state, pulse operator, and detection state use each nucleus's magnetogyric ratio relative to proton, matching the source's stated high-field pre-polarisation model.

## Inputs and units

- `spin_system` is the Spinach system structure, including spin magnetogyric ratios and spin count.
- `parameters.sweep` is the acquisition spectral-window width in hertz; its reciprocal sets the evolution timestep.
- `parameters.npoints` is a positive integer acquisition point count.
- `parameters.detection` must be `'uniaxial'` or `'quadrature'`.
- `parameters.flip_angle` is a real numeric scalar in radians, specified for protons. The weighted pulse operator scales the response for other nuclei by their gamma ratio relative to proton.
- `H`, `R`, and `K` are numeric matrices of equal dimensions: the Hamiltonian, relaxation superoperator, and kinetics superoperator supplied by the context function.

## Pulse and acquisition

The Liouvillian is assembled as `H + 1i*R + 1i*K`. With `weights = spin_system.inter.gammas/spin('1H')`, the routine forms the initial density operator as the weighted sum of per-spin `Lz` states, the detection coil as the weighted sum of per-spin `L+` states, and the pulse operator as the weighted sum of per-spin `Ly` operators. It applies the pulse with `step(spin_system,Sy,rho,parameters.flip_angle)`, then calls `evolution` with timestep `1/parameters.sweep`, interval count `parameters.npoints-1`, and mode `'observable'`.

For `'quadrature'`, the returned complex FID is left unchanged. For `'uniaxial'`, the routine takes its real part, removing imaginary-channel information used for frequency-sign discrimination.

## Guardrails

The routine checks that `H`, `R`, and `K` are numeric matrices with matching dimensions; sweep is a positive real scalar; point count is a positive integer; detection mode is one of the two listed strings; and flip angle is a real numeric scalar. These checks do not validate physical consistency of the supplied Spinach system or superoperators, and the real-scalar comparisons do not explicitly reject non-finite values.

## Links

- MATLAB source: https://github.com/IlyaKuprov/Spinach/blob/main/experiments/zulf/zerofield.m
- [Spinach Wiki: zerofield.m](https://spindynamics.org/wiki/index.php?title=zerofield.m)
