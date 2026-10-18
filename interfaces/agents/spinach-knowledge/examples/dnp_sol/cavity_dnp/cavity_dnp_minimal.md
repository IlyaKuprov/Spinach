# cavity_dnp_minimal

This example illustrates the preparation, selective pulses, and cavity reset in Fig. 1 and Appendix D of Lorriaux et al., arXiv:2610.19180v1. It is not an ensemble fit.

After adding Spinach and this directory to the MATLAB path, run:

```matlab
result=cavity_dnp_minimal();
```

The calculation requests 64 workers, subject to the parallel profile limit. The example writes no files and opens no figures.

## Model

The system contains an electron, 209Bi (I=9/2), and a truncated bosonic cavity. All twenty donor levels are retained in `sphten-liouv`. The magnetic field is 65.37 mT, the isotropic hyperfine coupling is 1.475 GHz, and the electron g factor is 2.00032; the nuclear Zeeman term uses Spinach's isotope table.

The cavity is resonant with the |6> to |15> transition. Its illustrative single-donor coupling is 100 Hz and its photon linewidth is 210 kHz. The Jaynes–Cummings interaction is selected with `inter.modes.strength.exchange='rwa'`; donor Zeeman frequencies remain in the laboratory-frame Hamiltonian. The photon bath is at zero temperature. Intrinsic spin relaxation, inhomogeneous broadening, collective enhancement, and resonator-frequency modulation are omitted.

## Preparation and pulse sequence

The donor starts at 25 mK. During preparation the cavity is disconnected, and `reduce` removes its spectator coordinate from the pulse calculation without discarding donor levels or changing formalism.

Repeated rectangular RF sweeps saturate all eight adjacent transitions in the lower donor manifold. After each pulse, `coherence` averages the donor RF phase. This sequence replaces the experimental band-limited-noise scrambling; the sweep residual and raw preparation trace are returned.

A selective microwave pulse shelves |6> population in |15>, then a selective RF pulse transfers |7> population to |6>. Both use `shaped_pulse_af` with the full donor drift and a linearly polarised electronic magnetic drive. The smaller direct nuclear drive is neglected. Donor phase averaging removes phase-dependent transients but retains electron–nuclear zero-quantum coherences. The calibrated rectangular pulses approximate, rather than reproduce, the experimental WURST waveforms; they are not ideal eigenstate swaps.

The cavity vacuum is attached using `BL1`, and the connected system is rebuilt. `evolution` follows spin-cavity exchange and photon loss for 20 s. Disconnecting the cavity during the pulses approximates experimental resonator detuning.

## Results and numerical settings

The returned trajectories are nuclear Iz, photon number, trace, and highest-retained-Fock-level population. `delta_iz` measures nuclear polarisation transfer from the saturated preparation, not a twofold enhancement of |6> population. The example does not use a coupled-eigenstate population label.

Traces use the native unit state, `state(spin_system,'E',1)`, and observables use left-product operators. The raw preparation trace is recorded before normalisation. Precision is set by `paranoia`; `spin_system.tols.small_matrix=size(H,1)+1` selects the matrix-exponential branch of Spinach's propagator for the bounded reset matrices.

The example uses three cavity levels, pulse Fokker–Planck rank 2, and 300 saturation sweeps. These settings were checked against cavity dimensions 4 and 5, rank 3, and 600 sweeps. The maximum nuclear Iz differences were below 3e-5, 1.4e-6, and 4.6e-7, respectively. The reset time grid is 0 to 20 s in 1 ms steps. These constants are visible near the top of the example rather than exposed as function arguments.

The same calculation was also checked with reset time-step refinements by factors of 2 and 5/4. Physical controls omitted the selective RF pulse or spin-cavity exchange, tested lossless conservation of donor projection plus photon number, resolved the photon transient, and compared one-photon decay with the exact exponential. Trace, observable bounds, saturation residual, and Fock-edge occupation were checked. No validation driver is shipped with this example.
