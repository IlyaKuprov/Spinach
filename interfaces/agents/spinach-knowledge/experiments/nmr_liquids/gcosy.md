# experiments/nmr_liquids/gcosy.m

- Signature: `fid=gcosy(spin_system,parameters,H,R,K)`

## Purpose

Horne-Morris gradient-selected COSY pulse sequence. Syntax: fid=gcosy(spin_system,parameters,H,R,K)

## Physical / mathematical content

- Liquid-state pulse sequence implementations. These are production experiment kernels that carry out coherence transfer, mixing, refocusing, decoupling, and indirect evolution on precomputed Hamiltonian/relaxation/kinetics operators.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.

## Parameters / inputs

- parameters.sweep sweep width in Hz
- parameters.npoints number of points for both dimensions
- parameters.spins nuclei on which the sequence runs,
- specified as {'1H'}, {'13C'}, etc.
- parameters.angle second pulse angle in radians, usu-
- ally pi/2, but also allows COSY45,
- COSY60, etc.
- parameters.g_amp gradient amplitude in Gauss/cm,
- defaults to 3
- parameters.g_dur gradient duration in seconds,
- defaults to 2e-3
- parameters.g_stab_del post-gradient stabilisation delay in
- seconds, defaults to 2e-4
- parameters.s_len active sample length in cm,
- defaults to 1.5
- parameters.pathway optional coherence pathway selection,
- either 'P', 'N', or 'P+N', defaults
- to 'P'
- H -Hamiltonian matrix, received from context function
- R -relaxation superoperator, received from context function
- K -kinetics superoperator, received from context function

## Outputs

- fid -two-dimensional free induction decay, or a structure
- with P-type fid.pos and N-type fid.neg fields in
- 'P+N' mode
- Note: the default P-type pathway uses opposite gradient signs
- and is less sensitive to mixing pulse phase errors. The
- N-type pathway uses equal gradient signs.
- Note: 'P+N' mode returns P-type and N-type components for
- echo/anti-echo recombination in phase-sensitive processing.

## Implementation structure

- Horne-Morris gradient-selected COSY pulse sequence. Syntax:
- fid=gcosy(spin_system,parameters,H,R,K)
- parameters.sweep sweep width in Hz
- parameters.npoints number of points for both dimensions
- parameters.spins nuclei on which the sequence runs,
- specified as {'1H'}, {'13C'}, etc.
- parameters.angle second pulse angle in radians, usu-
- ally pi/2, but also allows COSY45,
- COSY60, etc.
- parameters.g_amp gradient amplitude in Gauss/cm,
- defaults to 3
- parameters.g_dur gradient duration in seconds,
