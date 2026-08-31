# experiments/esr_dipolar/deer_3p_soft_deer.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/experiments/esr_dipolar/deer_3p_soft_deer.m`
- Signature: `echo_stack=deer_3p_soft_deer(spin_system,parameters,H,R,K)`
- Total lines: 236

## Purpose

Three-pulse DEER/PELDOR pulse sequence. The sequence uses soft pulses computed with the Fokker-Planck formalism. Syntax: echo_stack=deer_3p_soft_deer(spin_system,parameters,H,R,K)

## Physical / mathematical content

- Dipolar ESR experiment implementations. The pulse logic resolves dipolar couplings by echo modulation, with selective excitation and time-domain accumulation.
- The file uses a Fokker-Planck-style enlarged state space in which spatial or orientational coordinates are promoted to extra dimensions and coupled to spin dynamics through differential operators.
- Propagation is accelerated with a Krylov-subspace method, replacing direct matrix exponentiation by projection into a much smaller Arnoldi/Lanczos-type subspace.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.
- The spin physics includes through-space magnetic dipole-dipole coupling, a rank-2 anisotropic interaction with strong orientation dependence and characteristic secular/non-secular structure.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- A Krylov-subspace or Arnoldi construction is used to avoid forming or exponentiating very large dense propagators directly.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 77-78: Move into adjoint representation if needed; implemented by `[spin_system,parameters,H,R,K]=sim2liouv(spin_system,parameters,H,R,K)`.
- Lines 80-81: Check consistency; implemented by `grumble(spin_system,parameters,H,R,K)`.
- Lines 83-84: Compose Liouvillian; implemented by `L=H+1i*R+1i*K`.
- Lines 86-87: Electron pulse operators; implemented by `Ep=operator(spin_system,'L+',parameters.spins{1})`.
- Lines 90-92: Frequency offsets; implemented by `parameters.pulse_frq=-spin_system.inter.magnet*spin('E')/(2*pi)- parameters.pulse_frq-parameters.offset`.
- Lines 93-96: First pulse; implemented by `rho=shaped_pulse_af(spin_system,L,Ex,Ey,parameters.rho0,parameters.pulse_frq(1),parameters.pulse_pwr(1), parameters.pulse_dur(1),parameters.pulse_phi(1), parameters.puls…`.
- Lines 97-98: Evolution; implemented by `stepsize=parameters.p1_p3_gap/parameters.p2_nsteps`.
- Lines 101-104: Second pulse; implemented by `rho_stack=shaped_pulse_af(spin_system,L,Ex,Ey,rho_stack,parameters.pulse_frq(2),parameters.pulse_pwr(2), parameters.pulse_dur(2),parameters.pulse_phi(2), parameters.puls…`.
- Lines 105-106: Evolution; implemented by `rho_stack(:,end:-1:1)=evolution(spin_system,L,[],rho_stack(:,end:-1:1),stepsize,parameters.p2_nsteps,'refocus')`.
- Lines 108-111: Third pulse; implemented by `rho_stack=shaped_pulse_af(spin_system,L,Ex,Ey,rho_stack,parameters.pulse_frq(3),parameters.pulse_pwr(3), parameters.pulse_dur(3),parameters.pulse_phi(3), parameters.puls…`.
- Lines 112-114: Evolve to the edge of the echo window; implemented by `echo_location=parameters.p1_p3_gap+parameters.pulse_dur(1)/2+ parameters.pulse_dur(2)-parameters.echo_time/2`.
- Lines 117-118: Detect the echo; implemented by `stepsize=parameters.echo_time/parameters.echo_npts`.

### Key state/data transformations

- Lines 78: computes `[spin_system,parameters,H,R,K]` using `[spin_system,parameters,H,R,K]=sim2liouv(spin_system,parameters,H,R,K)`.
- Lines 84: computes `L` using `L=H+1i*R+1i*K`.
- Lines 87: computes `Ep` using `Ep=operator(spin_system,'L+',parameters.spins{1})`.
- Lines 88: computes `Ex` using `Ex=(Ep+Ep')/2; Ey=(Ep-Ep')/2i`.
- Lines 91-92: computes `parameters.pulse_frq` using `parameters.pulse_frq=-spin_system.inter.magnet*spin('E')/(2*pi)- parameters.pulse_frq-parameters.offset`.
- Lines 94-96: computes `rho` using `rho=shaped_pulse_af(spin_system,L,Ex,Ey,parameters.rho0,parameters.pulse_frq(1),parameters.pulse_pwr(1), parameters.pulse_dur(1),parameters.pulse_phi(1), parameters.puls…`.
- Lines 98: computes `stepsize` using `stepsize=parameters.p1_p3_gap/parameters.p2_nsteps`.
- Lines 99: computes `rho_stack` using `rho_stack=evolution(spin_system,L,[],rho,stepsize,parameters.p2_nsteps,'trajectory')`.
- Lines 106: computes `rho_stack(:,end:-1:1)` using `rho_stack(:,end:-1:1)=evolution(spin_system,L,[],rho_stack(:,end:-1:1),stepsize,parameters.p2_nsteps,'refocus')`.
- Lines 113-114: computes `echo_location` using `echo_location=parameters.p1_p3_gap+parameters.pulse_dur(1)/2+ parameters.pulse_dur(2)-parameters.echo_time/2`.
- Lines 119-120: computes `echo_stack` using `echo_stack=evolution(spin_system,L,parameters.coil,rho_stack, stepsize,parameters.echo_npts,'observable')`.

### Local helper functions

- Line 125: `grumble()` — `function grumble(spin_system,parameters,H,R,K)`.
  - Representative operation: `if (~isnumeric(H))||(~isnumeric(R))||(~isnumeric(K))|| (~ismatrix(H))||(~ismatrix(R))||(~ismatrix(K))`.
  - Representative operation: `(~ismatrix(H))||(~ismatrix(R))||(~ismatrix(K))`.

## Parameters / inputs

- parameters.pulse_frq -frequencies for the three
- pulses, Hz
- parameters.pulse_pwr -power levels for the three
- pulses, rad/s
- parameters.pulse_dur -durations for the three
- pulses, seconds
- parameters.pulse_phi -initial phases for the three
- pulses, radians
- parameters.pulse_rnk -Fokker-Planck ranks for the
- three pulses
- parameters.p1_p3_gap -time between the first and the
- third pulses, seconds
- parameters.p2_nsteps -number of second pulse posi-
- tions in the interval between
- the first and the third pulse
- parameters.echo_time -time to sample around the ex-
- pected echo position
- parameters.echo_npts -number of points in the echo
- discretization
- parameters.rho0 -initial state
- parameters.coil -detection state
- parameters.spins -irradiated spins, normally {'E'}
- parameters.method -soft puse propagation method,
- 'expv' for Krylov propagation,
- 'expm' for exponential propa-
- gation, 'evolution' for Spin-
- ach evolution function
- H -Hamiltonian matrix, received from context function
- R -relaxation superoperator, received from context function
- K -kinetics superoperator, received from context function

## Outputs

- echo_stack -DEER echo stack, a matrix with p2_nsteps echoes
- with echo_npts points each
- Note: for the method, start with 'expm', change to 'expv' if the
- calculation runs out of memory, and use 'evolution' as the
- last resort.
- Note: simulated echoes tend to be sharp and hard to catch becau-
- se simulation does not have distributions in experimental
- parameters. Fourier transforming the echo prior to integ-
- ration is recommended.
- Note: the time in the DEER trace refers to the second pulse inser-
- tion point, after end of the first pulse.

## Implementation structure

- Three-pulse DEER/PELDOR pulse sequence. The sequence uses soft
- pulses computed with the Fokker-Planck formalism. Syntax:
- echo_stack=deer_3p_soft_deer(spin_system,parameters,H,R,K)
- parameters.pulse_frq -frequencies for the three
- pulses, Hz
- parameters.pulse_pwr -power levels for the three
- pulses, rad/s
- parameters.pulse_dur -durations for the three
- pulses, seconds
- parameters.pulse_phi -initial phases for the three
- pulses, radians
- parameters.pulse_rnk -Fokker-Planck ranks for the

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `sim2liouv()`, `grumble()`, `operator()`, `spin()`, `shaped_pulse_af()`, `evolution()`, `rho_stack()`, `ismatrix()`, `all()`, `ismember()`, `isfield()`, `isscalar()`, `any()`, `iscell()`, `ischar()`.
