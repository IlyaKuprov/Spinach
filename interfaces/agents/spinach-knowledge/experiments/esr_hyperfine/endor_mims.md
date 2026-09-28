# experiments/esr_hyperfine/endor_mims.m

- Signature: `fid=endor_mims(spin_system,parameters,H,R,K)`

## Purpose

Mims ENDOR pulse sequence with ideal hard pulses. Syntax: fid=endor_mims(spin_system,parameters,H,R,K)

## Physical / mathematical content


- The Mims ENDOR sequence selects electron zero-order coherence, applies a nuclear pulse, and detects the resulting electron coherence; nuclear coherence orders are selected during the pulse pathway.
- The evolution generator combines the Hamiltonian, relaxation, and kinetics terms as `L = H + iR + iK`.

## Numerical / algorithmic content


- The indirect time-domain trajectory is sampled at spacing `1/parameters.sweep`; `parameters.npoints` controls the number of points.
- The routine propagates with `step()` and `evolution()` and returns the resulting time-domain ENDOR signal; FFT processing is downstream, not performed here.

## Parameters / inputs

- parameters.sweep nuclear frequency sweep width, Hz
- parameters.npoints number of fid points to be computed
- parameters.tau stimulated echo time, seconds
- H -Hamiltonian matrix, received from context function
- R -relaxation superoperator, received from context function
- K -kinetics superoperator, received from context function

## Outputs

- fid -free induction decay whose Fourier transform is the
- Mims ENDOR signal

## Implementation structure


- Converts the system to the adjoint representation when required, validates the Liouville-space inputs, and constructs electron and nuclear pulse operators.
- Applies the Mims pulse and coherence-selection sequence, evolves the indirect dimension, and returns the sampled trajectory.
