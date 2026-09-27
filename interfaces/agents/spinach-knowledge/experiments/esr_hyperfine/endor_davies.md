# experiments/esr_hyperfine/endor_davies.m

- Signature: `answer=endor_davies(spin_system,parameters,H,R,K)`

## Purpose


Simulates a Davies ENDOR sequence with explicit soft electron and nuclear pulses, including orientation-selection effects. The soft pulses use the Fokker-Planck formalism.

## Physical / mathematical content


- The sequence compares RF-on and RF-off branches; both undergo the electron-pulse sequence, while only the RF-on branch receives the nuclear pulse. Orientation selection is represented through soft pulses in the Fokker-Planck formalism.
- For each nuclear frequency, the detected signal is the RF-on amplitude divided by the RF-off reference amplitude.

## Numerical / algorithmic content


- Pulse evolution combines the Hamiltonian, relaxation, and kinetics terms as `L = H + iR + iK`; shaped pulses are propagated in the Fokker-Planck formalism.
- An optional electron spin echo is included when `parameters.tau` is nonzero, and the nuclear-frequency scan is evaluated in parallel.

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


- Converts to the adjoint representation when needed, validates the inputs, and constructs the electron and nuclear pulse operators.
- Applies the electron pulses, evaluates the RF-on and RF-off branches for each nuclear frequency, and returns their detected-amplitude ratio.
