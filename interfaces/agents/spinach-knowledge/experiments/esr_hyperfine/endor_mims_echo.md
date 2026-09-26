# experiments/esr_hyperfine/endor_mims_echo.m

- Signature: `stim_echo=endor_mims_echo(spin_system,parameters,H,R,K)`

## Purpose

Stimulated echo diagnostics for the Mims ENDOR sequence. Syntax: stim_echo=endor_mims_echo(spin_system,parameters,H,R,K)

## Physical / mathematical content

- Hyperfine ESR experiment implementations. These sequences probe coupled electron-nuclear dynamics through ENDOR or HYSCORE-type manipulations of coherence pathways.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.

## Parameters / inputs

- parameters.spins -working spins, normally {'E'}; spe-
- cify multiplicity if electron spin
- is not 1/2, for example {'7E'} for
- gadolinium
- parameters.electrons -a vector of integers specifying
- which spins in sys.isotopes are
- electrons
- parameters.tau -the delay between the first two
- 90-degree pulses of the Mims
- ENDOR sequence, seconds; 200e-9
- is typical
- parameters.n_dur -duration of the nuclear pulse,
- which is NOT ACTUALLY APPLIED
- within this pulse sequence,
- seconds; 50e-6 is typical
- parameters.nsteps -number of time steps to make
- in the detection period, which
- runs from 0 to 2*paramters.tau
- H -Hamiltonian matrix, received from the context
- function, normally powder() in this case
- R -relaxation superoperator, received from the context
- function, normally powder() in this case
- K -kinetics superoperator, received from the context
- function, normally powder() in this case

## Outputs

- stim_echo -stimulated echo seen in the Mims ENDOR sequ-
- ence in the absence of the nuclear RF pulse

## Implementation structure

- Stimulated echo diagnostics for the Mims ENDOR sequence. Syntax:
- stim_echo=endor_mims_echo(spin_system,parameters,H,R,K)
- parameters.spins -working spins, normally {'E'}; spe-
- cify multiplicity if electron spin
- is not 1/2, for example {'7E'} for
- gadolinium
- parameters.electrons -a vector of integers specifying
- which spins in sys.isotopes are
- electrons
- parameters.tau -the delay between the first two
- 90-degree pulses of the Mims
- ENDOR sequence, seconds; 200e-9
