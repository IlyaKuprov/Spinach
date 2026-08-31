# experiments/esr_dipolar/deer_4p_soft_deer.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/experiments/esr_dipolar/deer_4p_soft_deer.m`
- Signature: `echo_stack=deer_4p_soft_deer(spin_system,parameters,H,R,K)`
- Total lines: 268

## Purpose

Four-pulse DEER/PELDOR pulse sequence. The sequence uses soft pulses computed with the Fokker-Planck formalism. Syntax: echo_stack=deer_4p_soft_deer(spin_system,parameters,H,R,K)

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

- Lines 83-84: Move into adjoint representation if needed; implemented by `[spin_system,parameters,H,R,K]=sim2liouv(spin_system,parameters,H,R,K)`.
- Lines 86-87: Check consistency; implemented by `grumble(spin_system,parameters,H,R,K)`.
- Lines 89-90: Compose Liouvillian; implemented by `L=H+1i*R+1i*K`.
- Lines 92-93: Electron pulse operators; implemented by `Ep=operator(spin_system,'L+',parameters.spins{1})`.
- Lines 96-98: Frequency offsets; implemented by `parameters.pulse_frq=-spin_system.inter.magnet*spin('E')/(2*pi)- parameters.pulse_frq-parameters.offset`.
- Lines 100-103: First pulse; implemented by `rho=shaped_pulse_af(spin_system,L,Ex,Ey,parameters.rho0,parameters.pulse_frq(1),parameters.pulse_pwr(1), parameters.pulse_dur(1),parameters.pulse_phi(1), parameters.puls…`.
- Lines 104-105: Evolution; implemented by `rho=evolution(spin_system,L,[],rho,parameters.p1_p2_gap,1,'final')`.
- Lines 107-110: Second pulse; implemented by `rho=shaped_pulse_af(spin_system,L,Ex,Ey,rho,parameters.pulse_frq(2),parameters.pulse_pwr(2), parameters.pulse_dur(2),parameters.pulse_phi(2), parameters.pulse_rnk(2),par…`.
- Lines 111-112: Evolution; implemented by `rho=evolution(spin_system,L,[],rho,parameters.p1_p2_gap+parameters.pulse_dur(1)/2,1,'final')`.
- Lines 114-115: Evolution; implemented by `stepsize=(parameters.p2_p4_gap-parameters.p1_p2_gap)/parameters.p3_nsteps`.
- Lines 118-121: Third pulse; implemented by `rho=shaped_pulse_af(spin_system,L,Ex,Ey,rho,parameters.pulse_frq(3),parameters.pulse_pwr(3), parameters.pulse_dur(3),parameters.pulse_phi(3), parameters.pulse_rnk(3),par…`.
- Lines 123-124: Evolution; implemented by `rho(:,end:-1:1)=evolution(spin_system,L,[],rho(:,end:-1:1),stepsize,parameters.p3_nsteps,'refocus')`.
- Lines 126-129: Fourth pulse; implemented by `rho=shaped_pulse_af(spin_system,L,Ex,Ey,rho,parameters.pulse_frq(4),parameters.pulse_pwr(4), parameters.pulse_dur(4),parameters.pulse_phi(4), parameters.pulse_rnk(4),par…`.
- Lines 130-132: Evolve to the start of the echo; implemented by `echo_location=(parameters.p2_p4_gap-parameters.p1_p2_gap)+ parameters.pulse_dur(3)-parameters.echo_time/2`.
- Lines 135-136: Sample the echo; implemented by `stepsize=parameters.echo_time/parameters.echo_npts`.

### Key state/data transformations

- Lines 84: computes `[spin_system,parameters,H,R,K]` using `[spin_system,parameters,H,R,K]=sim2liouv(spin_system,parameters,H,R,K)`.
- Lines 90: computes `L` using `L=H+1i*R+1i*K`.
- Lines 93: computes `Ep` using `Ep=operator(spin_system,'L+',parameters.spins{1})`.
- Lines 94: computes `Ex` using `Ex=(Ep+Ep')/2; Ey=(Ep-Ep')/2i`.
- Lines 97-98: computes `parameters.pulse_frq` using `parameters.pulse_frq=-spin_system.inter.magnet*spin('E')/(2*pi)- parameters.pulse_frq-parameters.offset`.
- Lines 101-103: computes `rho` using `rho=shaped_pulse_af(spin_system,L,Ex,Ey,parameters.rho0,parameters.pulse_frq(1),parameters.pulse_pwr(1), parameters.pulse_dur(1),parameters.pulse_phi(1), parameters.puls…`.
- Lines 115: computes `stepsize` using `stepsize=(parameters.p2_p4_gap-parameters.p1_p2_gap)/parameters.p3_nsteps`.
- Lines 124: computes `rho(:,end:-1:1)` using `rho(:,end:-1:1)=evolution(spin_system,L,[],rho(:,end:-1:1),stepsize,parameters.p3_nsteps,'refocus')`.
- Lines 131-132: computes `echo_location` using `echo_location=(parameters.p2_p4_gap-parameters.p1_p2_gap)+ parameters.pulse_dur(3)-parameters.echo_time/2`.
- Lines 137-138: computes `echo_stack` using `echo_stack=evolution(spin_system,L,parameters.coil,rho, stepsize,parameters.echo_npts,'observable')`.

### Local helper functions

- Line 143: `grumble()` — `function grumble(spin_system,parameters,H,R,K)`.
  - Representative operation: `if (~isnumeric(H))||(~isnumeric(R))||(~isnumeric(K))|| (~ismatrix(H))||(~ismatrix(R))||(~ismatrix(K))`.
  - Representative operation: `(~ismatrix(H))||(~ismatrix(R))||(~ismatrix(K))`.

## Parameters / inputs

- parameters.pulse_frq -frequencies for the four
- pulses, Hz
- parameters.pulse_pwr -power levels for the four
- pulses, rad/s
- parameters.pulse_dur -durations for the four
- pulses, seconds
- parameters.pulse_phi -initial phases for the four
- pulses, radians
- parameters.pulse_rnk -Fokker-Planck ranks for the
- four pulses
- parameters.p1_p2_gap -time between the end of the
- first and the start of the
- second pulse, seconds
- parameters.p2_p4_gap -time between the end of the
- second the start of the third
- pulse, seconds
- parameters.p3_nsteps -number of third pulse posi-
- tions in the interval between
- the first echo and the fourth
- pulse
- parameters.echo_time -time to sample around the ex-
- pected second echo position
- parameters.echo_npts -number of points in the second
- echo discretization
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

- echo_stack -DEER echo stack, a matrix with p3_nsteps echoes
- with echo_npts points each
- Note: for the method, start with 'expm', change to 'expv' if the
- calculation runs out of memory, and use 'evolution' as the
- last resort.
- Note: simulated echoes tend to be sharp and hard to catch becau-
- se simulation does not have distributions in experimental
- parameters. Fourier transforming the echo prior to integ-
- ration is recommended.
- Note: the time in the DEER trace refers to the third pulse inser-
- tion point, after end of the second pulse.

## Implementation structure

- Four-pulse DEER/PELDOR pulse sequence. The sequence uses soft
- pulses computed with the Fokker-Planck formalism. Syntax:
- echo_stack=deer_4p_soft_deer(spin_system,parameters,H,R,K)
- parameters.pulse_frq -frequencies for the four
- pulses, Hz
- parameters.pulse_pwr -power levels for the four
- pulses, rad/s
- parameters.pulse_dur -durations for the four
- pulses, seconds
- parameters.pulse_phi -initial phases for the four
- pulses, radians
- parameters.pulse_rnk -Fokker-Planck ranks for the

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `sim2liouv()`, `grumble()`, `operator()`, `spin()`, `shaped_pulse_af()`, `evolution()`, `rho()`, `ismatrix()`, `all()`, `ismember()`, `isfield()`, `isscalar()`, `any()`, `iscell()`, `ischar()`.
