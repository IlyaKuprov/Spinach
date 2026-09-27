# experiments/esr_hyperfine/endor_mims_echo.m

- Signature: `stim_echo=endor_mims_echo(spin_system,parameters,H,R,K)`

## Purpose

Stimulated echo diagnostics for the Mims ENDOR sequence. Syntax: stim_echo=endor_mims_echo(spin_system,parameters,H,R,K)

## Physical / mathematical content


- Implements an electron stimulated-echo diagnostic for Mims ENDOR, using electron Lz as the initial state and electron L+ for detection. It does not apply a nuclear RF pulse; the nuclear-pulse-duration parameter is used as an evolution interval.
- The sequence uses ideal electron pi/2 rotations about x, x, and y, with the specified delays between them.

## Numerical / algorithmic content


- The state is propagated through the pulse and delay sequence with the Liouvillian `L = H + iR + iK`, and the echo is sampled over the configured time points.
- Propagation uses `step()` and `evolution()`; the returned signal is the detected time-domain echo, not a HYSCORE dataset.

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


- Converts to Liouville/adjoint representation as needed, checks dimensions and parameters, and constructs electron pulse and detection operators.
- Applies the electron rotations and free-evolution intervals, then returns the sampled stimulated echo.
