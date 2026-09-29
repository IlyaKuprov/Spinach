# experiments/esr_hyperfine/endor_davies.m

- MATLAB implementation: [experiments/esr_hyperfine/endor_davies.m](https://github.com/IlyaKuprov/Spinach/blob/main/experiments/esr_hyperfine/endor_davies.m)

Signature: answer=endor_davies(spin_system,parameters,H,R,K)

## Purpose and physical sequence

This is a Davies ENDOR simulation with explicit shaped electron and nuclear pulses. The soft-pulse treatment uses the Fokker–Planck formalism and can represent orientation selection. The code first applies an electron pi pulse. For each nuclear RF frequency it compares an RF-on nuclear pulse with a same-duration zero-power reference, applies an electron pi/2 pulse to each branch, and returns the ratio of the two coil-detected signals. If parameters.tau is supplied and nonzero, both branches also pass through a spin-echo stage: free evolution for tau, an electron pi pulse, and a second tau period. Thus the output is a simulated ENDOR response over the entries of n_frq, not a magnetic-field axis or measured spectrum.

## Inputs and required settings

- spin_system; H, R, and K are the spin system and context-provided Hamiltonian, relaxation, and kinetics matrices. The matrices must have matching dimensions; the routine moves to the adjoint representation if needed.
- Electron pulse: parameters.e_frq in Hz, e_phi in radians, e_pwr in rad/s, e_dur in seconds, and positive-integer Fokker–Planck cutoff e_rnk. The source specifies e_dur as the electron pi-pulse duration and obtains the pi/2 pulse by using half that duration.
- Nuclear pulse: parameters.n_frq is a vector of RF frequencies in Hz; n_phi is in radians, n_pwr in rad/s, n_dur in seconds, and n_rnk is the positive-integer Fokker–Planck cutoff.
- parameters.method selects the method passed to shaped_pulse_af. parameters.spins lists the irradiated spins with electron first and nucleus second. parameters.offset supplies electron and nuclear transmitter offsets in Hz.
- parameters.rho0 is the initial state and parameters.coil is the detection state. Optional parameters.tau is a non-negative spin-echo delay in seconds; omitting it skips the echo stage.

## Output axis and numerical cautions

answer has the same shape as parameters.n_frq; each entry is the coil signal with nuclear RF divided by the zero-power reference signal at that frequency. No example numeric RF frequency or pulse duration is given in the source, so none is invented here. Increase the electron and nuclear Fokker–Planck ranks and the spherical-grid size until the output converges, as the source notes. The function is an ENDOR spin-dynamics simulation; it does not calculate DNP/hyperpolarisation, image encoding, or a field sweep.

Source: https://spindynamics.org/wiki/index.php?title=endor_davies.m
