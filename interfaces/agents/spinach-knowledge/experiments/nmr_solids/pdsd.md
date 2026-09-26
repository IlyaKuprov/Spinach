# experiments/nmr_solids/pdsd.m

- Signature: `fid=pdsd(spin_system,parameters,H,R,K)`

## Purpose

A simplified model of the PDSD experiment using NOESY type quadrature detection and phase cycle. To be cal- led from the singlerot context. Syntax: fid=pdsd(spin_system,parameters,H,R,K)

## Physical / mathematical content

- Solid-state pulse sequence implementations. The core ingredients are anisotropic Hamiltonians, rotor synchronisation, cross-polarisation, recoupling/decoupling, and powder or rotor-stack propagation.
- Orientation or trajectory averaging is performed numerically, so grid design, weights, and integration error control matter directly to accuracy and runtime.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- Numerical integration over angles or geometry is part of the implementation, so point placement and weights are as important as the local Hamiltonian calculations.

## Parameters / inputs

- spin_system -Spinach spin system object
- parameters.sweep -sweep width in Hz
- parameters.npoints -two-element vector giving the
- number of complex points in the
- indirect and direct dimensions
- parameters.tmix -mixing time in seconds
- parameters.rate -MAS rate in Hz, used to set
- proton irradiation power
- parameters.spc_dim -spatial dimension of the MAS
- problem, received from the
- context function
- H, R, K -Hamiltonian, relaxation, and
- kinetics superoperators, recei-
- ved from the context function

## Outputs

- fid.cos, fid.sin -States quadrature components
- of the 2D PDSD spectrum

## Implementation structure

- A simplified model of the PDSD experiment using NOESY
- type quadrature detection and phase cycle. To be cal-
- led from the singlerot context. Syntax:
- fid=pdsd(spin_system,parameters,H,R,K)
- spin_system -Spinach spin system object
- parameters.sweep -sweep width in Hz
- parameters.npoints -two-element vector giving the
- number of complex points in the
- indirect and direct dimensions
- parameters.tmix -mixing time in seconds
- parameters.rate -MAS rate in Hz, used to set
- proton irradiation power
