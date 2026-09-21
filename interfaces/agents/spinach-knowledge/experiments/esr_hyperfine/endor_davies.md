# experiments/esr_hyperfine/endor_davies.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/experiments/esr_hyperfine/endor_davies.m`
- Signature: `answer=endor_davies(spin_system,parameters,H,R,K)`
- Total lines: 269

## Purpose

Davies ENDOR sequence with explicit soft pulses and all of the atten- dant effects, such as orientation selection. Soft pulses are simula- ted using the Fokker-Planck formalism. Syntax: answer=endor_davies(spin_system,parameters,H,R,K)

## Physical / mathematical content

- Hyperfine ESR experiment implementations. These sequences probe coupled electron-nuclear dynamics through ENDOR or HYSCORE-type manipulations of coherence pathways.
- The file uses a Fokker-Planck-style enlarged state space in which spatial or orientational coordinates are promoted to extra dimensions and coupled to spin dynamics through differential operators.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Parameters / inputs

- The following parameters refer to the electron pi pulse. The duration
- of the electron pi/2 pulse is obtained by halving parameters.e_dur:
- parameters.e_frq -frequency of the electron pulse, Hz
- parameters.e_phi -phase of the electron pulse, rad
- parameters.e_pwr -power of the electron pulse, rad/s
- parameters.e_dur -duration of the electron pulse, s
- parameters.e_rnk -Fokker-Planck cut-off rank for
- the electron pulse
- The following parameters refer to the nuclei pulse:
- parameters.n_frq -vector of frequencies for the nuclei
- pulse, in Hz. The answer is returned
- as a vector of the same dimension.
- parameters.n_phi -phase of the nuclei pulse, rad
- parameters.n_pwr -power of the nuclei pulse, rad/s
- parameters.n_dur -duration of the nuclei pulse, s
- parameters.n_rnk -Fokker-Planck cut-off rank for
- the nuclei pulse
- parameters.method -method to use during the call
- to shaped_pulse_af()
- parameters.spins -irradiated spins, electron first,
- nucleus second
- parameters.offset -transmitter offsets for the electron
- and the nucleus pulses, Hz
- parameters.rho0 -initial state
- parameters.coil -detection state
- parameters.tau -optional spin echo delay, seconds
- H -Hamiltonian matrix, received from context function
- R -relaxation superoperator, received from context function
- K -kinetics superoperator, received from context function

## Outputs

- answer -amplitude detected on the coil state for each
- frequency of the nuclear pulse
- Note: Fokker-Planck ranks should be increased until convergence is
- achieved in the output. The same applies to the size of the
- spherical grid.

## Implementation structure

- Davies ENDOR sequence with explicit soft pulses and all of the atten-
- dant effects, such as orientation selection. Soft pulses are simula-
- ted using the Fokker-Planck formalism. Syntax:
- answer=endor_davies(spin_system,parameters,H,R,K)
- The following parameters refer to the electron pi pulse. The duration
- of the electron pi/2 pulse is obtained by halving parameters.e_dur:
- parameters.e_frq -frequency of the electron pulse, Hz
- parameters.e_phi -phase of the electron pulse, rad
- parameters.e_pwr -power of the electron pulse, rad/s
- parameters.e_dur -duration of the electron pulse, s
- parameters.e_rnk -Fokker-Planck cut-off rank for
- the electron pulse

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `sim2liouv()`, `grumble()`, `operator()`, `spin()`, `shaped_pulse_af()`, `isfield()`, `evolution()`, `answer()`, `ismatrix()`, `all()`, `ismember()`, `ischar()`, `iscell()`, `cellfun()`, `isrow()`, `isscalar()`.
