# experiments/nmr_liquids/hetcor.m

- Signature: `fid=hetcor(spin_system,parameters,H,R,K)`

## Purpose

Magnitude-mode HETCOR pulse sequence from:

## Physical / mathematical content

- Liquid-state pulse sequence implementations. These are production experiment kernels that carry out coherence transfer, mixing, refocusing, decoupling, and indirect evolution on precomputed Hamiltonian/relaxation/kinetics operators.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.

## Syntax

```matlab
fid=hetcor(spin_system,parameters,H,R,K)
```

## Parameters / inputs

- parameters.sweep [F1 F2] sweep widths in the
- two frequency directions, Hz
- parameters.npoints [F1 F2] numbers of points in
- the two time directions in fid
- parameters.spins {F1 F2} nuclei (e.g. {'1H','13C'})
- parameters.decouple list of nuclei that detection
- time decoupling should be applied
- to -cell array of strings, e.g.
- {'1H','15N'})
- parameters.J working scalar coupling, Hz
- H -Hamiltonian matrix, received from context function
- R -relaxation superoperator, received from context function
- K -kinetics superoperator, received from context function

## Outputs

- fid -two-dimensional free induction decay for magnitude-
- mode processing
- Note: natural abundance experiments should make use of the iso-
- tope dilution functionality. See dilute.m function.
- Note: the fixed transfer delays are delta_2=1/(2J) and
- delta_3=1/(3J), using the absolute value of J.

## Implementation structure

- Magnitude-mode HETCOR pulse sequence from:
- fid=hetcor(spin_system,parameters,H,R,K)
- parameters.sweep [F1 F2] sweep widths in the
- two frequency directions, Hz
- parameters.npoints [F1 F2] numbers of points in
- the two time directions in fid
- parameters.spins {F1 F2} nuclei (e.g. {'1H','13C'})
- parameters.decouple list of nuclei that detection
- time decoupling should be applied
- to -cell array of strings, e.g.
- {'1H','15N'})
- parameters.J working scalar coupling, Hz
