# experiments/esr_hyperfine/endor_mims_ideal.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/experiments/esr_hyperfine/endor_mims_ideal.m`
- Signature: `endor_spec=endor_mims_ideal(spin_system,parameters,H,R,K)`
- Total lines: 196

## Purpose

Mims ENDOR sequence with ideal electron pulses. Syntax: endor_spec=endor_mims_ideal(spin_system,parameters,H,R,K)

## Physical / mathematical content

- Hyperfine ESR experiment implementations. These sequences probe coupled electron-nuclear dynamics through ENDOR or HYSCORE-type manipulations of coherence pathways.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

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
- seconds; 50e-6 is typical
- parameters.n_frq -nuclear pulse frequency offsets
- parameters.rf_b1_field -RF B1 field strength
- parameters.n_rnk -nuclear pulse grid rank
- parameters.nuclei -a vector of integers specifying
- which spins in sys.isotopes are
- nuclei to irradiate
- H -Hamiltonian matrix, received from the context
- function, normally powder() in this case
- R -relaxation superoperator, received from the context
- function, normally powder() in this case
- K -kinetics superoperator, received from the context
- function, normally powder() in this case

## Outputs

- endor_spec -Mims ENDOR spectrum, a vector of the same
- size as parameters.n_frq

## Implementation structure

- Mims ENDOR sequence with ideal electron pulses. Syntax:
- endor_spec=endor_mims_ideal(spin_system,parameters,H,R,K)
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

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `sim2liouv()`, `grumble()`, `operator()`, `state()`, `nut_freqs()`, `step()`, `rho()`, `shaped_pulse_af()`, `ismatrix()`, `all()`, `ismember()`, `isfield()`, `iscell()`, `ischar()`, `isscalar()`, `isrow()`.
