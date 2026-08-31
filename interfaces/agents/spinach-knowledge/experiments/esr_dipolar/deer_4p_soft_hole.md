# experiments/esr_dipolar/deer_4p_soft_hole.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/experiments/esr_dipolar/deer_4p_soft_hole.m`
- Signature: `fids=deer_4p_soft_hole(spin_system,parameters,H,R,K)`
- Total lines: 200

## Purpose

Pulse diagnostics for the four-pulse DEER/PELDOR pulse sequen- ce. This function shows how soft pulses affect the magnetisati- on of the sample. It is a hypothetical experiment where a soft pulse specified by the user is performed, immediately followed by an ideal pi/2 pulse on all spins followed by infinite-band- width time-domain detection. Syntax: fids=deer_4p_soft_hole(spin_system,parameters,H,R,K)

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

- Lines 69-70: Move into adjoint representation if needed; implemented by `[spin_system,parameters,H,R,K]=sim2liouv(spin_system,parameters,H,R,K)`.
- Lines 72-73: Check consistency; implemented by `grumble(spin_system,parameters,H,R,K)`.
- Lines 75-76: Compose Liouvillian; implemented by `L=H+1i*R+1i*K`.
- Lines 78-79: Pulse operators; implemented by `Ep=operator(spin_system,'L+',parameters.spins{1})`.
- Lines 82-84: Frequency offsets; implemented by `parameters.pulse_frq=-spin_system.inter.magnet*spin('E')/(2*pi)- parameters.pulse_frq-parameters.offset`.
- Lines 86-89: Soft pulses; implemented by `rho1=shaped_pulse_af(spin_system,L,Ex,Ey,parameters.rho0,parameters.pulse_frq(1),parameters.pulse_pwr(1), parameters.pulse_dur(1),parameters.pulse_phi(1), parameters.pul…`.
- Lines 99-100: Hard pulse; implemented by `parameters.rho0=step(spin_system,Ey,[rho1 rho2 rho3 rho4],pi/2)`.
- Lines 102-103: Acquisition; implemented by `fids=acquire(spin_system,parameters,H,R,K)`.

### Key state/data transformations

- Lines 70: computes `[spin_system,parameters,H,R,K]` using `[spin_system,parameters,H,R,K]=sim2liouv(spin_system,parameters,H,R,K)`.
- Lines 76: computes `L` using `L=H+1i*R+1i*K`.
- Lines 79: computes `Ep` using `Ep=operator(spin_system,'L+',parameters.spins{1})`.
- Lines 80: computes `Ex` using `Ex=(Ep+Ep')/2; Ey=(Ep-Ep')/2i`.
- Lines 83-84: computes `parameters.pulse_frq` using `parameters.pulse_frq=-spin_system.inter.magnet*spin('E')/(2*pi)- parameters.pulse_frq-parameters.offset`.
- Lines 87-89: computes `rho1` using `rho1=shaped_pulse_af(spin_system,L,Ex,Ey,parameters.rho0,parameters.pulse_frq(1),parameters.pulse_pwr(1), parameters.pulse_dur(1),parameters.pulse_phi(1), parameters.pul…`.
- Lines 90-92: computes `rho2` using `rho2=shaped_pulse_af(spin_system,L,Ex,Ey,parameters.rho0,parameters.pulse_frq(2),parameters.pulse_pwr(2), parameters.pulse_dur(2),parameters.pulse_phi(2), parameters.pul…`.
- Lines 93-95: computes `rho3` using `rho3=shaped_pulse_af(spin_system,L,Ex,Ey,parameters.rho0,parameters.pulse_frq(3),parameters.pulse_pwr(3), parameters.pulse_dur(3),parameters.pulse_phi(3), parameters.pul…`.
- Lines 96-98: computes `rho4` using `rho4=shaped_pulse_af(spin_system,L,Ex,Ey,parameters.rho0,parameters.pulse_frq(4),parameters.pulse_pwr(4), parameters.pulse_dur(4),parameters.pulse_phi(4), parameters.pul…`.
- Lines 100: computes `parameters.rho0` using `parameters.rho0=step(spin_system,Ey,[rho1 rho2 rho3 rho4],pi/2)`.
- Lines 103: computes `fids` using `fids=acquire(spin_system,parameters,H,R,K)`.

### Local helper functions

- Line 108: `grumble()` — `function grumble(spin_system,parameters,H,R,K)`.
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
- parameters.offset -receiver offset for the time
- domain detection, Hz
- parameters.sweep -sweep width for time domain
- detection, Hz
- parameters.npoints -number of points in the free
- induction decay
- parameters.spins -irradiated spins, normally {'E'}
- parameters.rho0 -initial state
- parameters.coil -detection state
- parameters.method -soft puse propagation method,
- 'expv' for Krylov propagation,
- 'expm' for exponential propa-
- gation, 'evolution' for Spin-
- ach evolution function
- H -Hamiltonian matrix, received from context function
- R -relaxation superoperator, received from context function
- K -kinetics superoperator, received from context function

## Outputs

- fids -four free induction decays that should be apo-
- dised and Fourier transformed
- Note: for the method, start with 'expm', change to 'expv' if the
- calculation runs out of memory, and use 'evolution' as the
- last resort.

## Implementation structure

- Pulse diagnostics for the four-pulse DEER/PELDOR pulse sequen-
- ce. This function shows how soft pulses affect the magnetisati-
- on of the sample. It is a hypothetical experiment where a soft
- pulse specified by the user is performed, immediately followed
- by an ideal pi/2 pulse on all spins followed by infinite-band-
- width time-domain detection. Syntax:
- fids=deer_4p_soft_hole(spin_system,parameters,H,R,K)
- parameters.pulse_frq -frequencies for the four
- pulses, Hz
- parameters.pulse_pwr -power levels for the four
- pulses, rad/s
- parameters.pulse_dur -durations for the four

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `sim2liouv()`, `grumble()`, `operator()`, `spin()`, `shaped_pulse_af()`, `step()`, `acquire()`, `ismatrix()`, `all()`, `ismember()`, `isfield()`, `isscalar()`, `any()`, `iscell()`, `ischar()`.
