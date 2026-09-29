# experiments/esr_hyperfine/hyscore.m

- MATLAB implementation: [experiments/esr_hyperfine/hyscore.m](https://github.com/IlyaKuprov/Spinach/blob/main/experiments/esr_hyperfine/hyscore.m)

Source: https://spindynamics.org/wiki/index.php?title=hyscore.m

Signature: fid=hyscore(spin_system,parameters,H,R,K).

This simulates a HYSCORE (hyperfine sublevel correlation) electron-spin-echo modulation experiment. It probes electron-nuclear couplings through two time dimensions of the echo modulation; the output is a free-induction decay, not a measured spectrum. The source attributes the implementation to Szosenfogel and Goldfarb, DOI https://doi.org/10.1080/00268979809483260.

The routine forms the Liouvillian L=H+1i*R+1i*K and electron Lx from the E raising/lowering operators. It applies an ideal x-axis pi/2 pulse, evolves for tau, applies another pi/2 pulse, and filters to zero electron coherence ({'E',0}) to select the relevant nuclear-modulation pathways. It records the indirect evolution trajectory, applies the third ideal electron pi pulse, propagates the detection state backward through tau, applies the corresponding backward pi/2 rotation, and records the direct-dimension observable trajectory.

Required inputs:

- parameters.nsteps: two positive integer counts, [n1 n2], for the indirect and direct time dimensions.
- parameters.sweep: scalar sweep width in Hz. The source uses 1/sweep as the time increment in each dimension.
- parameters.tau: echo delay in seconds.
- parameters.rho0: initial state vector.
- parameters.coil: detection state vector.
- H, R, and K: same-size matrices defining Hamiltonian, relaxation, and kinetics contributions; the routine checks their matrix dimensions against one another and the states.

Return value: fid is the two-dimensional time-domain free-induction decay indexed by the two nsteps dimensions, with increments 1/sweep seconds. Fourier-transforming this array gives the HYSCORE spectrum; the transform is not performed by this routine.

Limit: all three electron pulses are ideal hard rotations. The source explicitly recommends replacing them with shaped_pulse_af() when soft pulses are required. No numerical defaults for tau, sweep, or nsteps are specified in the source.
