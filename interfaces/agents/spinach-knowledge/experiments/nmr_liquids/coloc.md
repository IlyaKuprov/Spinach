# experiments/nmr_liquids/coloc.m

- Signature: `fid=coloc(spin_system,parameters,H,R,K)`

## Purpose

COLOC NMR pulse sequence from: Implemented as shown in Fig 1b, without the dashed pulses during the delta(2) period. Delta(1) defaults to half of the maximum F1 evolution time implied by the sweep width, delta(2) must be spe- cified. Syntax: fid=coloc(spin_system,parameters,H,R,K)

## Physical / mathematical content

- Liquid-state pulse sequence implementations. These are production experiment kernels that carry out coherence transfer, mixing, refocusing, decoupling, and indirect evolution on precomputed Hamiltonian/relaxation/kinetics operators.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.

## Parameters / inputs

- parameters.sweep [F1 F2] sweep widths, Hz
- parameters.npoints [F1 F2] numbers of points
- parameters.spins {F1 F2} nuclei (e.g. '13C','1H')
- parameters.delta2 COLOC delta2 (see the paper),
- typically 40e-3 seconds
- parameters.delta1 optional COLOC delta1 delay, seconds;
- must be at least half of maximum t1
- H -Hamiltonian matrix, received from context function
- R -relaxation superoperator, received from context function
- K -kinetics superoperator, received from context function

## Outputs

- fid -free induction decay for magnitude mode processing
- Note: natural abundance simulations should make use of the isotope
- dilution functionality. See dilute.m function.

## Implementation structure

- COLOC NMR pulse sequence from:
- Implemented as shown in Fig 1b, without the dashed pulses during
- the delta(2) period. Delta(1) defaults to half of the maximum F1
- evolution time implied by the sweep width, delta(2) must be spe-
- cified. Syntax:
- fid=coloc(spin_system,parameters,H,R,K)
- parameters.sweep [F1 F2] sweep widths, Hz
- parameters.npoints [F1 F2] numbers of points
- parameters.spins {F1 F2} nuclei (e.g. '13C','1H')
- parameters.delta2 COLOC delta2 (see the paper),
- typically 40e-3 seconds
- parameters.delta1 optional COLOC delta1 delay, seconds;
