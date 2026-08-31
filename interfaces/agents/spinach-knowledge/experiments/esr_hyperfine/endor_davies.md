# experiments/esr_hyperfine/endor_davies.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/experiments/esr_hyperfine/endor_davies.m`
- Signature: `answer=endor_davies(spin_system,parameters,H,R,K)`
- Total lines: 269

## Purpose

Davies ENDOR sequence with explicit soft pulses and all of the atten- dant effects, such as orientation selection. Soft pulses are simula- ted using the Fokker-Planck formalism. Syntax: answer=endor_davies(spin_system,parameters,H,R,K)

## Physical / mathematical content

- Hyperfine ESR experiment implementations. These sequences probe coupled electron-nuclear dynamics through ENDOR or HYSCORE-type manipulations of coherence pathways.
- The file uses a Fokker-Planck-style enlarged state space in which spatial or orientational coordinates are promoted to extra dimensions and coupled to spin dynamics through differential operators.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 74-75: Move into adjoint representation if needed; implemented by `[spin_system,parameters,H,R,K]=sim2liouv(spin_system,parameters,H,R,K)`.
- Lines 77-78: Check consistency; implemented by `grumble(spin_system,parameters,H,R,K)`.
- Lines 80-81: Compose Liouvillian; implemented by `L=H+1i*R+1i*K`.
- Lines 83-84: Pulse operators; implemented by `Ep=operator(spin_system,'L+',parameters.spins{1})`.
- Lines 87-88: Pulse operators; implemented by `Np=operator(spin_system,'L+',parameters.spins{2})`.
- Lines 91-92: Soft pulse frequencies; implemented by `parameters.e_frq=parameters.e_frq-parameters.offset(1)`.
- Lines 95-99: Pi pulse on the electron; implemented by `rho0=shaped_pulse_af(spin_system,L,Ex,Ey,parameters.rho0, parameters.e_frq, parameters.e_pwr,parameters.e_dur, parameters.e_phi,parameters.e_rnk, parameters.method)`.
- Lines 101-102: Preallocate the answer; implemented by `answer=zeros(size(parameters.n_frq),'like',1i)`.
- Lines 104-105: Loop over radiofrequency offsets; implemented by `parfor n=1:numel(parameters.n_frq)`.
- Lines 107-110: Either pi pulse on the nuclei or a delay of the same length for reference; implemented by `rho_a=shaped_pulse_af(spin_system,L,Nx,Ny,rho0,parameters.n_frq(n),parameters.n_pwr, parameters.n_dur, parameters.n_phi, parameters.n_rnk, parameters.method)`.
- Lines 115-118: Pi/2 pulse on the electron; implemented by `rho_a=shaped_pulse_af(spin_system,L,Ex,Ey,rho_a,parameters.e_frq, parameters.e_pwr, parameters.e_dur/2,parameters.e_phi, parameters.e_rnk, parameters.method)`.
- Lines 123-124: Optionally, run the spin echo stage; implemented by `if isfield(parameters,'tau')&&(parameters.tau~=0)`.
- Lines 126-127: First evolution period; implemented by `rho_a=evolution(spin_system,L,[],rho_a,parameters.tau,1,'final')`.
- Lines 130-133: Pi pulse on the electron; implemented by `rho_a=shaped_pulse_af(spin_system,L,Ex,Ey,rho_a,parameters.e_frq,parameters.e_pwr, parameters.e_dur,parameters.e_phi, parameters.e_rnk,parameters.method)`.
- Lines 138-139: Second evolution period; implemented by `rho_a=evolution(spin_system,L,[],rho_a,parameters.tau,1,'final')`.
- Lines 144-145: Return the ratio of RF-on and RF-off; implemented by `answer(n)=(parameters.coil'*rho_a)/(parameters.coil'*rho_b)`.

### Control flow inferred from the code

- Line 105: `parfor` loop over `n=1:numel(parameters.n_frq)`.
- Line 124: conditional branch on `isfield(parameters,'tau')&&(parameters.tau~=0)`.

### Key state/data transformations

- Lines 75: computes `[spin_system,parameters,H,R,K]` using `[spin_system,parameters,H,R,K]=sim2liouv(spin_system,parameters,H,R,K)`.
- Lines 81: computes `L` using `L=H+1i*R+1i*K`.
- Lines 84: computes `Ep` using `Ep=operator(spin_system,'L+',parameters.spins{1})`.
- Lines 85: computes `Ex` using `Ex=(Ep+Ep')/2; Ey=(Ep-Ep')/2i`.
- Lines 88: computes `Np` using `Np=operator(spin_system,'L+',parameters.spins{2})`.
- Lines 89: computes `Nx` using `Nx=(Np+Np')/2; Ny=(Np-Np')/2i`.
- Lines 92: computes `parameters.e_frq` using `parameters.e_frq=parameters.e_frq-parameters.offset(1)`.
- Lines 93: computes `parameters.n_frq` using `parameters.n_frq=parameters.n_frq+spin(parameters.spins{2})*spin_system.inter.magnet/(2*pi)`.
- Lines 96-99: computes `rho0` using `rho0=shaped_pulse_af(spin_system,L,Ex,Ey,parameters.rho0, parameters.e_frq, parameters.e_pwr,parameters.e_dur, parameters.e_phi,parameters.e_rnk, parameters.method)`.
- Lines 102: computes `answer` using `answer=zeros(size(parameters.n_frq),'like',1i)`.
- Lines 108-110: computes `rho_a` using `rho_a=shaped_pulse_af(spin_system,L,Nx,Ny,rho0,parameters.n_frq(n),parameters.n_pwr, parameters.n_dur, parameters.n_phi, parameters.n_rnk, parameters.method)`.
- Lines 111-113: computes `rho_b` using `rho_b=shaped_pulse_af(spin_system,L,Nx,Ny,rho0,parameters.n_frq(n),0*parameters.n_pwr, parameters.n_dur, parameters.n_phi, parameters.n_rnk, parameters.method)`.
- Lines 145: computes `answer(n)` using `answer(n)=(parameters.coil'*rho_a)/(parameters.coil'*rho_b)`.

### Local helper functions

- Line 152: `grumble()` — `function grumble(spin_system,parameters,H,R,K)`.
  - Representative operation: `if (~isnumeric(H))||(~isnumeric(R))||(~isnumeric(K))|| (~ismatrix(H))||(~ismatrix(R))||(~ismatrix(K))`.
  - Representative operation: `(~ismatrix(H))||(~ismatrix(R))||(~ismatrix(K))`.

## Parameters / inputs

- The following parameters refer to the electron pi pulse. The duration
- of the electron pi/2 pulse is obtained by halving parameters.e_dur:
- parameters.e_frq -frequency of the electron pulse, Hz
- parameters.e_phi -phase of the electron pulse, rad
- parameters.e_pwr -power of the electron pulse, rad/s
- parameters.e_dur -duration of the electron pulse, s
- parameters.e_rnk -Fokker-Planck cut-off rank for
- the electron pulse
- The following parameters refer to the nuclei pulse:
- parameters.n_frq -vector of frequencies for the nuclei
- pulse, in Hz. The answer is returned
- as a vector of the same dimension.
- parameters.n_phi -phase of the nuclei pulse, rad
- parameters.n_pwr -power of the nuclei pulse, rad/s
- parameters.n_dur -duration of the nuclei pulse, s
- parameters.n_rnk -Fokker-Planck cut-off rank for
- the nuclei pulse
- parameters.method -method to use during the call
- to shaped_pulse_af()
- parameters.spins -irradiated spins, electron first,
- nucleus second
- parameters.offset -transmitter offsets for the electron
- and the nucleus pulses, Hz
- parameters.rho0 -initial state
- parameters.coil -detection state
- parameters.tau -optional spin echo delay, seconds
- H -Hamiltonian matrix, received from context function
- R -relaxation superoperator, received from context function
- K -kinetics superoperator, received from context function

## Outputs

- answer -amplitude detected on the coil state for each
- frequency of the nuclear pulse
- Note: Fokker-Planck ranks should be increased until convergence is
- achieved in the output. The same applies to the size of the
- spherical grid.

## Implementation structure

- Davies ENDOR sequence with explicit soft pulses and all of the atten-
- dant effects, such as orientation selection. Soft pulses are simula-
- ted using the Fokker-Planck formalism. Syntax:
- answer=endor_davies(spin_system,parameters,H,R,K)
- The following parameters refer to the electron pi pulse. The duration
- of the electron pi/2 pulse is obtained by halving parameters.e_dur:
- parameters.e_frq -frequency of the electron pulse, Hz
- parameters.e_phi -phase of the electron pulse, rad
- parameters.e_pwr -power of the electron pulse, rad/s
- parameters.e_dur -duration of the electron pulse, s
- parameters.e_rnk -Fokker-Planck cut-off rank for
- the electron pulse

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `sim2liouv()`, `grumble()`, `operator()`, `spin()`, `shaped_pulse_af()`, `isfield()`, `evolution()`, `answer()`, `ismatrix()`, `all()`, `ismember()`, `ischar()`, `iscell()`, `cellfun()`, `isrow()`, `isscalar()`.
