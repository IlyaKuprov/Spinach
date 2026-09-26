# experiments/bruker/hmqcetgp.m

- Signature: `fid=hmqcetgp(spin_system,parameters,H,R,K)`

## Purpose

Echo/antiecho gradient-selected HMQC pulse sequence, based on the Bruker hmqcetgp pulse program and the standard HMQC sequence from: The gradient selection is represented analytically by coherence order selection statements

## Physical / mathematical content

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.

## Syntax

```matlab
fid=hmqcetgp(spin_system,parameters,H,R,K)
```

## Parameters / inputs

- parameters.sweep [F1 F2] sweep widths, Hz
- parameters.npoints [F1 F2] numbers of points
- parameters.spins {F1 F2} nuclei (e.g. '13C','1H')
- parameters.decouple_f2 nuclei to decouple in F2, e.g.
- {'15N','13C'}
- parameters.decouple_f1 nuclei that receive midpoint
- 180-degree refocusing pulses in
- F1, e.g. {'1H'}
- parameters.J working scalar coupling, Hz
- H -Hamiltonian matrix, received from context function
- R -relaxation superoperator, received from context function
- K -kinetics superoperator, received from context function

## Outputs

- fid.pos,fid.neg -echo and antiecho components of the
- signal.
- Note: natural abundance simulations should make use of the isotope
- dilution functionality. See dilute.m function.

## Implementation structure

- Echo/antiecho gradient-selected HMQC pulse sequence, based on the
- Bruker hmqcetgp pulse program and the standard HMQC sequence from:
- The gradient selection is represented analytically by coherence order
- selection statements
- fid=hmqcetgp(spin_system,parameters,H,R,K)
- parameters.sweep [F1 F2] sweep widths, Hz
- parameters.npoints [F1 F2] numbers of points
- parameters.spins {F1 F2} nuclei (e.g. '13C','1H')
- parameters.decouple_f2 nuclei to decouple in F2, e.g.
- {'15N','13C'}
- parameters.decouple_f1 nuclei that receive midpoint
- 180-degree refocusing pulses in
