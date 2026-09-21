# kernel/eigenfields/eigenfields.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/eigenfields/eigenfields.m`
- Signature: `tran=eigenfields(spin_system,parameters,Hz,Hc,Hmw)`
- Total lines: 514

## Purpose

Computes resonance fields. For a Hamiltonian Hc+B*Hz, returns all magnetic fields B for which the difference between two eigenvalues of Hc+B*Hz is equal to the frequency provided, and the transition moment across the specified operator Hmw is significant. Syntax: tran=eigenfields(spin_system,parameters,Hz,Hc,Hmw)

## Physical / mathematical content

- Eigenfield utilities. These files analyse field-dependent eigenstructure and resonance conditions, linking Hamiltonian spectra to magnetic-field sweeps and transition behaviour.

## Numerical / algorithmic content

- An eigenvalue problem is solved or analysed, so the file is extracting spectra, stationary states, avoided crossings, or modal structure from the effective Hamiltonian or superoperator.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `size()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Parameters / inputs

- Hz -field-dependent part of the laboratory frame Hamil-
- tonian operator (Hilbert space) or commutation supe-
- roperator (Liouville space), normalised to 1 Tesla
- Hc -field-independent part of the laboratory frame Hamil-
- tonian operator (Hilbert space) or commutation supe-
- roperator (Liouville space), containing couplings and
- offsets
- Hmw -observable operator (Hilbert space) or observable
- vector (Liouville space), without the amplitude
- prefactor
- parameters.window -magnet field window, Tesla
- parameters.mw_freq -microwave frequency, Hz
- parameters.orientation -three Euler angles in radians
- specifying the system orientation
- parameters.tm_tol -relative transition moment
- tolerance
- parameters.pp_tol -peak position tolerance in Tesla,
- this should be much smaller than
- the typical line width
- parameters.fwhm -transition full width at half
- maximum, Tesla
- parameters.rspt_order -perturbation theory order to use
- to account for the off-diagonal
- part of the Hamiltonian, Inf for
- exact diagonalisation

## Outputs

- tran.tf -vector of transition fields in Tesla
- tran.tm -vector of transition moments
- tran.tw -vector of transition FWHMs in Tesla
- tran.pd -vector of energy level population differences
- tran.ti -transition identity array, one row per transition
- tran.tj -vector of scaled field-sweep Jacobians

## Implementation structure

- Computes resonance fields. For a Hamiltonian Hc+B*Hz, returns all
- magnetic fields B for which the difference between two eigenvalues
- of Hc+B*Hz is equal to the frequency provided, and the transition
- moment across the specified operator Hmw is significant. Syntax:
- tran=eigenfields(spin_system,parameters,Hz,Hc,Hmw)
- Hz - field-dependent part of the laboratory frame Hamil-
- tonian operator (Hilbert space) or commutation supe-
- roperator (Liouville space), normalised to 1 Tesla
- Hc - field-independent part of the laboratory frame Hamil-
- roperator (Liouville space), containing couplings and
- offsets
- Hmw - observable operator (Hilbert space) or observable

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `spin()`, `rspt_eig()`, `false()`, `any()`, `converged()`, `new_conv()`, `true()`, `herm_spline()`, `new_grid()`, `new_dE()`, `new_LP()`, `new_E()`, `new_T()`, `new_V()`, `state_ovlp()`.
