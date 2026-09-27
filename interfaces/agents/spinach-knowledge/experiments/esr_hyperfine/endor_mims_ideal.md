# experiments/esr_hyperfine/endor_mims_ideal.m

- Signature: `endor_spec=endor_mims_ideal(spin_system,parameters,H,R,K)`

## Purpose

Mims ENDOR sequence with ideal electron pulses. Syntax: endor_spec=endor_mims_ideal(spin_system,parameters,H,R,K)

## Physical / mathematical content


- Implements Mims ENDOR: electron polarization is stored along `Lz`, followed by electron pi/2 pulses, a nuclear RF pulse, and electron-coherence detection through `L+`.
- Nuclear pulse operators are assembled from x- and y-axis rotations weighted by the gyromagnetic ratios of the selected nuclei; the evolution generator is `L = H + iR + iK`.

## Numerical / algorithmic content


- The sequence propagates the density operator through the configured pulse and delay intervals and evaluates the detected response across the nuclear-frequency points.
- The implementation uses `step()` and `evolution()` for propagation; it returns the ENDOR response and does not perform an FFT or other spectrum post-processing.

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


- Converts to the adjoint representation when needed and checks dimensions, Liouville formalism, and parameter shapes before constructing pulse operators.
- Builds electron and gyromagnetic-ratio-weighted nuclear operators, applies the pulse-delay sequence, and returns the detected response.
