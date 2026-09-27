# experiments/nmr_solids/cn2d_sq.m

- Signature: `fid=cn2d_sq(spin_system,parameters,H,R,K)`

## Purpose

Single-quantum version of the 13C-detected 14N-13C MAS 2D correlation experiment described by Jarvis, Haies, Williamson, and Carravetta ([paper](http://dx.doi.org/10.1039/c3cp50787d)).

## Physical / mathematical content

`parameters.spins` specifies 14N first and 13C second. The implementation builds pulse operators for those spins and applies coherence-order selection for the single-quantum sequence.

## Numerical / algorithmic content

The source evolves separate cosine and sine quadrature branches, using finite-duration 14N RF pulses and trajectory evolution over the sweep increments, with 13C pulses and refocusing in the sequence.

## Parameters / inputs

- parameters.spins -isotopes to which the sequence is
- applied, specified as a cell array
- with 14N first, and 13C second
- parameters.spc_dim -Fokker-Planck spatial dimension
- parameters.sweep -sweep widths in the two dimensions, Hz
- parameters.npoints -numbers of points in the two dimensions
- parameters.rho0 -initial state
- parameters.coil -detection state
- parameters.rf_pwr -RF power on 14N, Hz
- parameters.rf_dur -RF pulse duration on 14N, seconds

## Outputs
- fid.sin
- fid.cos -sine and cosine components
- of the States quadrature
