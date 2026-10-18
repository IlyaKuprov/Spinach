# Full low-field Bi:Si cavity-assisted DNP cycle

Run `result=cavity_dnp_minimal(3)` and `validation=validate_cavity_dnp(result)` after adding Spinach and this directory to the MATLAB path. The example requests 64 workers; the configured parallel profile may cap that request. Neither function writes files or opens figures.

## Physical model

The particles are the actual electron, **209Bi (I=9/2)**, and a truncated bosonic cavity. All twenty donor levels are retained. The field is 65.37 mT, the isotropic hyperfine coupling is 1.475 GHz, and the electron g factor is 2.00032; the nuclear gyromagnetic ratio is Spinach's native isotope value. The laboratory-frame Zeeman terms and complete isotropic hyperfine interaction remain present during preparation, pulses, and reset. The cavity frequency is computed from the native donor spectrum for the paper's |6>↔|15> transition.

This single-donor time-domain illustration follows Fig. 1 and the reset logic of Appendix D in Lorriaux et al., arXiv:2610.19180v1. It is not an ensemble fit. The bare electron–cavity exchange coupling is an illustrative 100 Hz, and the photon linewidth is 210 kHz. Native `inter.modes.strength.exchange='rwa'` selects the Jaynes–Cummings branch; it does not remove donor Zeeman frequencies. The photon bath is approximated as zero-temperature. Intrinsic spin relaxation, inhomogeneous broadening, collective enhancement, and resonator-frequency modulation are omitted.

## Physical preparation and pulses

1. Native `equilibrium` prepares the full donor at 25 mK. While the cavity is disconnected, its identity coordinate is only a computational spectator. Native `reduce` identifies the exact invariant donor subspace; this is a projection within `sphten-liouv`, not a change of formalism or a truncation of donor levels.
2. Finite rectangular RF pulses address all eight adjacent transitions in the lower donor manifold. Native `coherence(...,{{[1 2],0}})` implements the uniform donor RF-phase ensemble analytically. Repeated sweeps saturate the manifold, replacing the experimental random-phase WURST scrambling with a simpler physical saturation sequence. The returned sweep residual and preparation trace are measured, and doubling the number of sweeps is tested.
3. A finite selective microwave pulse shelves |6> population into |15>, followed by a finite RF pulse transferring |7> population into |6>. Both use native `shaped_pulse_af` with the full donor drift and a linearly polarized electronic magnetic drive. The same native total-donor coherence-order selection suppresses phase-dependent transients while preserving electron–nuclear zero-quantum coherences. Durations are calibrated rectangular-pulse approximations, not ideal eigenstate swaps or experimental WURST waveforms. The smaller direct nuclear magnetic-drive term is neglected.
4. Native bosonic `BL1` syntax attaches vacuum; the connected Hamiltonian and relaxation are rebuilt from the interaction specification. Native multichannel `evolution` then follows 20 s of laboratory-frame spin–cavity exchange and photon loss. The cavity is disconnected during the pulses, approximating the experimental resonator detuning.

Only eigenvalues of the native **Liouville-space left-product Hamiltonian** are used to determine the donor transition frequencies. No eigenvectors, coupled-eigenstate projectors, density matrices, Kronecker products, custom dissipators, or custom propagators are constructed. Native pulse maps are reused for saturation, and native `coherence` performs the analytical phase cycling; no averaged propagator or phase-cycle implementation is constructed.

## Observables and normalization

The outputs are nuclear Iz, intracavity photon occupation, trace, and highest-retained-Fock-level occupation. The reported `delta_iz` is the change from the physically saturated donor preparation, not a claim of twofold |6> population enhancement. Spinach has no coupled-eigenstate label used by this example.

`state(s,'E',1)` supplies one identity trace functional; `state(s,'E','all')` would sum an identity for every particle and is not used. Native left-product observable duals, `operator(...,'left')'*unit`, give physical traces without assuming that an arbitrary Euclidean state overlap is an expectation value. The preparation's raw trace is recorded before its final normalization. Native `paranoia` precision settings are selected. For the bounded reset matrices, `s.tols.small_matrix=size(H,1)+1` selects the MATLAB matrix-exponential branch already implemented inside Spinach's `propagator`; the example constructs no propagator itself.

## Validation

The validator compares Fock cutoffs 3, 4, and 5, increases the pulse Fokker–Planck rank, doubles the saturation length, and tests both half-size and non-dyadic (0.8 ms) reset timesteps. The non-dyadic comparison does not merely reuse the same binary scaling-and-squaring path. It independently removes the RF pulse or the exchange coupling, checks lossless conservation of total spin projection plus photon number, resolves the actual donor-generated photon transient, and verifies native one-photon decay against exp(-kappa*t). It checks trace, physical observable bounds, saturation residual, and Fock-edge occupation. The full native laboratory Hamiltonian, not an effective three-level rate equation, drives every reset trajectory.
