# experiments/esr_dipolar/deer_3p_soft_diag.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/experiments/esr_dipolar/deer_3p_soft_diag.m`
- Signature: `deer_3p_soft_diag(spin_system,parameters)`
- Total lines: 272

## Purpose

Complete set of simulations related to three-pulse DEER. Runs pulse diagnostics, which is followed by echo diagnostics, which is follow- ed by DEER simulation. Syntax: deer_3p_soft_diag(spin_system,parameters)

## Physical / mathematical content

- Dipolar ESR experiment implementations. The pulse logic resolves dipolar couplings by echo modulation, with selective excitation and time-domain accumulation.
- The file uses a Fokker-Planck-style enlarged state space in which spatial or orientational coordinates are promoted to extra dimensions and coupled to spin dynamics through differential operators.
- Propagation is accelerated with a Krylov-subspace method, replacing direct matrix exponentiation by projection into a much smaller Arnoldi/Lanczos-type subspace.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.
- The spin physics includes through-space magnetic dipole-dipole coupling, a rank-2 anisotropic interaction with strong orientation dependence and characteristic secular/non-secular structure.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.
- A Krylov-subspace or Arnoldi construction is used to avoid forming or exponentiating very large dense propagators directly.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 88-89: Check consistency; implemented by `grumble(parameters)`.
- Lines 91-92: Run pulse diagnostics experiment; implemented by `fids=powder(spin_system,@deer_3p_soft_hole,parameters,parameters.assumptions)`.
- Lines 94-95: Apodisation and Fourier transform; implemented by `fid_a=apodisation(spin_system,fids(:,1),{{'crisp'}}); spectrum_a=fftshift(fft(fid_a,parameters.zerofill))`.
- Lines 100-101: Plotting; implemented by `kfigure(); scale_figure([2.0 1.5])`.
- Lines 108-109: Get DEER echo stack; implemented by `echo_stack=powder(spin_system,@deer_3p_soft_deer,parameters,parameters.assumptions)`.
- Lines 111-113: Time axis for the echo window; implemented by `echo_axis=1e9*linspace(-parameters.echo_time/2, parameters.echo_time/2,parameters.echo_npts+1)`.
- Lines 115-116: Time axis for the DEER trace; implemented by `deer_axis=1e6*linspace(0,parameters.p1_p3_gap,parameters.p2_nsteps+1)`.
- Lines 118-119: Axis set for the echo stack; implemented by `[deer_axis_2d,echo_axis_2d]=meshgrid(deer_axis,echo_axis)`.
- Lines 121-122: Plot the echo stack; implemented by `kfigure(); surf(deer_axis_2d,echo_axis_2d,real(echo_stack))`.
- Lines 126-127: Extract and phase the echo modulation; implemented by `[deer_echoes,deer_sigmas,deer_traces]=svd(echo_stack)`.
- Lines 131-132: Plot echo components; implemented by `kfigure(); plot(echo_axis,real(deer_echoes(:,1:3)))`.
- Lines 140-141: Plot DEER components; implemented by `kfigure(); plot(deer_axis,real(deer_traces(:,1:3)))`.

### Key state/data transformations

- Lines 92: computes `fids` using `fids=powder(spin_system,@deer_3p_soft_hole,parameters,parameters.assumptions)`.
- Lines 95: computes `fid_a` using `fid_a=apodisation(spin_system,fids(:,1),{{'crisp'}}); spectrum_a=fftshift(fft(fid_a,parameters.zerofill))`.
- Lines 96: computes `fid_b` using `fid_b=apodisation(spin_system,fids(:,2),{{'crisp'}}); spectrum_b=fftshift(fft(fid_b,parameters.zerofill))`.
- Lines 97: computes `fid_c` using `fid_c=apodisation(spin_system,fids(:,3),{{'crisp'}}); spectrum_c=fftshift(fft(fid_c,parameters.zerofill))`.
- Lines 98: computes `fid_d` using `fid_d=apodisation(spin_system,fids(:,4),{{'crisp'}}); spectrum_d=fftshift(fft(fid_d,parameters.zerofill))`.
- Lines 109: computes `echo_stack` using `echo_stack=powder(spin_system,@deer_3p_soft_deer,parameters,parameters.assumptions)`.
- Lines 112-113: computes `echo_axis` using `echo_axis=1e9*linspace(-parameters.echo_time/2, parameters.echo_time/2,parameters.echo_npts+1)`.
- Lines 116: computes `deer_axis` using `deer_axis=1e6*linspace(0,parameters.p1_p3_gap,parameters.p2_nsteps+1)`.
- Lines 119: computes `[deer_axis_2d,echo_axis_2d]` using `[deer_axis_2d,echo_axis_2d]=meshgrid(deer_axis,echo_axis)`.
- Lines 127: computes `[deer_echoes,deer_sigmas,deer_traces]` using `[deer_echoes,deer_sigmas,deer_traces]=svd(echo_stack)`.
- Lines 128: computes `deer_echoes` using `deer_echoes=deer_echoes*deer_sigmas/deer_traces(1)`.
- Lines 129: computes `deer_traces` using `deer_traces=deer_traces*deer_sigmas/deer_traces(1)`.

### Local helper functions

- Line 152: `grumble()` — `function grumble(parameters)`.
  - Representative operation: `if ~isfield(parameters,'pulse_frq')`.
  - Representative operation: `error('pulse frequencies must be specified in parameters.pulse_frq field.')`.

## Parameters / inputs

- parameters.pulse_frq -frequencies for the three
- pulses, Hz
- parameters.pulse_pwr -power levels for the three
- pulses, Hz
- parameters.pulse_dur -durations for the three
- pulses, seconds
- parameters.pulse_phi -initial phases for the three
- pulses, radians
- parameters.pulse_rnk -Fokker-Planck ranks for the
- three pulses
- parameters.offset -receiver offset for the time
- domain detection, Hz
- parameters.sweep -sweep width for time domain
- detection, Hz
- parameters.npoints -number of points in the free
- induction decay
- parameters.zerofill -length of the zero-filled FFT
- parameters.rho0 -initial state
- parameters.coil -detection state
- parameters.p1_p3_gap -time between the first and the
- third pulses, seconds
- parameters.p2_nsteps -number of second pulse posi-
- tions in the interval between
- the first and the third pulse
- parameters.echo_time -time to sample around the ex-
- pected echo position
- parameters.echo_npts -number of points in the echo
- discretization
- parameters.method -soft puse propagation method,
- 'expv' for Krylov propagation,
- 'expm' for exponential propa-
- gation, 'evolution' for Spin-
- ach evolution function
- parameters.assumptions -Hamiltonian generation assump-
- tions, use 'deer' to keep two-
- electron flip-flop terms and
- 'deer-zz' to drop them

## Outputs

- Figure 1: pulse diagnostics
- Figure 2: DEER echo stack
- Figure 3: principal components of the stack, echo
- Figure 4: principal components of the stack, DEER
- Note: for the method, start with 'expm', change to 'expv' if the
- calculation runs out of memory, and use 'evolution' as the
- last resort.
- Note: simulated echoes tend to be sharp and hard to catch becau-
- se simulation does not have distributions in experimental
- parameters. Fourier transforming the echo prior to integ-
- ration is recommended.
- Note: the time in the DEER trace refers to the second pulse inser-
- tion point, after end of first pulse.

## Implementation structure

- Complete set of simulations related to three-pulse DEER. Runs pulse
- diagnostics, which is followed by echo diagnostics, which is follow-
- ed by DEER simulation. Syntax:
- deer_3p_soft_diag(spin_system,parameters)
- parameters.pulse_frq -frequencies for the three
- pulses, Hz
- parameters.pulse_pwr -power levels for the three
- parameters.pulse_dur -durations for the three
- pulses, seconds
- parameters.pulse_phi -initial phases for the three
- pulses, radians
- parameters.pulse_rnk -Fokker-Planck ranks for the

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `powder()`, `apodisation()`, `fids()`, `fftshift()`, `kfigure()`, `scale_figure()`, `subplot()`, `plot_1d()`, `ktitle()`, `kylabel()`, `kxlabel()`, `deer_traces()`, `deer_echoes()`, `klegend()`, `isfield()`.
