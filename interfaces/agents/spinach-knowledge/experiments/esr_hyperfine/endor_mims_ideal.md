# experiments/esr_hyperfine/endor_mims_ideal.m

- MATLAB implementation: [experiments/esr_hyperfine/endor_mims_ideal.m](https://github.com/IlyaKuprov/Spinach/blob/main/experiments/esr_hyperfine/endor_mims_ideal.m)

Source: https://spindynamics.org/wiki/index.php?title=endor_mims_ideal.m

Signature: endor_spec=endor_mims_ideal(spin_system,parameters,H,R,K).

This is a simulated Mims ENDOR experiment for electron-nuclear hyperfine-coupled systems. Starting from electron longitudinal polarisation Lz, the routine applies ideal electron x-axis pi/2 (90°) pulses, separated by the Mims delay tau. At each entry of parameters.n_frq, it applies a shaped nuclear RF pulse to the selected nuclei twice, with RF on and RF off, and subtracts the off response as background. It then applies an electron y-axis minus-pi/2 pulse, evolves for a second tau, and detects the electron coherence with L+. The propagated Liouvillian is assembled from the supplied Hamiltonian, relaxation, and kinetics matrices as H + 1i*R + 1i*K.

Required inputs:

- parameters.spins: working-spin labels; the source gives {'E'} as the usual spin-1/2 electron example and {'7E'} for a non-spin-1/2 electron such as gadolinium.
- parameters.electrons: indices of electron spins in sys.isotopes; these define the electron pulse and detection operators.
- parameters.nuclei: indices of nuclei in sys.isotopes to irradiate.
- parameters.tau: Mims interpulse delay in seconds. The source identifies 200e-9 s (200 ns) as typical.
- parameters.n_dur: nuclear pulse duration in seconds; the source identifies 50e-6 s (50 microseconds) as typical.
- parameters.n_frq: array of nuclear RF frequency offsets in Hz, which is the sampled axis of the output.
- parameters.rf_b1_field: RF B1 field strength in tesla; the source uses it with nuclear gyromagnetic ratios to form nutation frequencies and specifies no default value.
- parameters.n_rnk: rank of the nuclear pulse grid (dimensionless numerical setting).
- H, R, and K: same-size matrices supplied by the calling context, normally powder() according to the source comments. They encode the Hamiltonian, relaxation, and kinetics contributions.

Return value: endor_spec is one signal value per entry of n_frq, in the same order. It is a simulated, RF-on-minus-RF-off complex coherence signal; the function does not claim a measured spectrum.

Limits: the electron pulses are ideal instantaneous rotations, whereas only the nuclear RF pulse is shaped. Interpretation therefore depends on the supplied spin system and context matrices, the chosen powder/orientation treatment, and the frequency grid; this routine itself does not add experimental data or perform a separate Fourier transform.
