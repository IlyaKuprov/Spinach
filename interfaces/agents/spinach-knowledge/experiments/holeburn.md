# experiments/holeburn.m

- Signature: fid=holeburn(spin_system,parameters,H,R,K)

## Purpose and physical scope

Models a spectral hole-burning sequence: a frequency-selective soft pulse acts on the supplied initial state, followed by a hard pi/2 observation pulse and free-induction-decay acquisition. The soft pulse is propagated in a rank-truncated Fokker–Planck representation. Hyperfine interactions can influence the result through the supplied Spinach Hamiltonian, but this function itself is not an ESEEM or ENDOR sequence and does not create a hyperfine model.

## Inputs and parameters

H, R and K are the Hamiltonian, relaxation and kinetics matrices supplied by the experiment context. The Liouville core accepts `sphten-liouv` and `zeeman-liouv` directly. A `zeeman-hilb` density-matrix context is also supported: `sim2liouv` converts the generators, basis and state-like `rho0`/`coil` fields to `zeeman-liouv` before the grumbler runs. Its combined generator is H+1i*R+1i*K.

- parameters.spins: the function uses the first spin label to build the pulse operator; the source code uses this field even though its header does not document it.
- parameters.pulse_frq: soft-pulse frequency in Hz.
- parameters.offset: receiver offset in Hz; the code subtracts it from pulse_frq before soft-pulse propagation.
- parameters.pulse_phi: soft-pulse phase in radians.
- parameters.pulse_pwr: soft-pulse amplitude in radians per second.
- parameters.pulse_dur: soft-pulse duration in seconds.
- parameters.pulse_rnk: Fokker–Planck cutoff rank. The source advises increasing the rank until the output converges.
- parameters.method: soft-pulse propagation choice: expv, expm, or evolution.
- parameters.rho0 and parameters.coil: initial state and detection state.
- parameters.sweep: acquisition sweep width in Hz; parameters.npoints is the number of FID points.

## Sequence and output

The routine forms the selected-spin raising operator, projects it into the enlarged space, and derives the x and y pulse operators. It applies shaped_pulse_af with the offset-corrected frequency, power, duration, phase, rank, and selected propagation method. It then applies a hard pi/2 step about the y operator and calls acquire with the resulting state. The only returned quantity is fid: the FID sampled at the requested point count and sweep width. No time-axis output is returned separately; the acquisition dwell is set by the sweep width.

The source provides pi/2 as the observation-pulse angle and recommends rank convergence, but gives no numeric pulse-frequency, power, duration, rank, or acquisition example. Those settings should be selected for the spin system and experiment rather than inferred here.

Source: https://github.com/IlyaKuprov/Spinach/blob/main/experiments/holeburn.m
