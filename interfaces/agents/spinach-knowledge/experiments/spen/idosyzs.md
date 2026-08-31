# experiments/spen/idosyzs.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/experiments/spen/idosyzs.m`
- Signature: `inten=idosyzs(spin_system,parameters,H,R,K,G,F)`
- Total lines: 206

## Purpose

A simplified model sequence of the ZS iDOSY pulse sequence. Syntax: inten=idosyzs(spin_system,parameters,H,R,K,G,F) This sequence must be called from the imaging() context, which would provide H, R, K, G and F. Parameters: parameters.rho0 -initial state parameters.coil -detection state parameters.spins -nuclei on which the sequence runs parameters.g_amp -gradient amplitude for diffusion encoding (T/m) parameters.sel_

## Physical / mathematical content

- SPEN experiment implementations. These files combine shaped pulses, gradients, spatial encoding, and often diffusion-aware propagation.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 51-52: Check consistency; implemented by `grumble(spin_system,parameters,H,R,K,G,F)`.
- Lines 54-55: Compose Liouvillian; implemented by `L=H+F+1i*R+1i*K`.
- Lines 57-58: Get pulse operators; implemented by `Lp=operator(spin_system,'L+',parameters.spins{1})`.
- Lines 62-63: Diffusion time evolution delay; implemented by `parameters.delta=parameters.delta_big-parameters.delta_sml-parameters.rf_dur`.
- Lines 65-66: Shaped pulse setup; implemented by `[amps,phases,~,~,scaling_factor]=read_wave(parameters.filename,parameters.pulse_npoints)`.
- Lines 69-70: Gradient amplitudes for the gradient during the selective pulse; implemented by `gradient_amplitudes=parameters.sel_g_amp*ones(parameters.pulse_npoints,1)`.
- Lines 72-74: Calculate the maximum RF field strength and calibrate the amplitude of the soft pulse (this calibration is the same as in BRUKER TopSpin); implemented by `gamma_B1=parameters.rf_phi/parameters.rf_dur`.
- Lines 78-79: Pulse transformation from (amplitude, phase) to (X,Y); implemented by `[Cx,Cy]=polar2cartesian(amps,phases)`.
- Lines 81-82: Apply the first pulse; implemented by `rho=step(spin_system,Ly,parameters.rho0,pi/2)`.
- Lines 84-85: Select coherence -1; implemented by `rho=coherence(spin_system,rho,{{parameters.spins{1},-1}})`.
- Lines 87-88: Evolve under the first positive gradient; implemented by `rho=evolution(spin_system,L+parameters.g_amp*G{1},[],rho,parameters.delta_sml,1,'final')`.
- Lines 90-91: Run the diffusion time evolution; implemented by `rho=evolution(spin_system,L,[],rho,parameters.delta,1,'final')`.
- Lines 93-95: Execution of the 180 degree shaped pulse; implemented by `rho=shaped_pulse_xy(spin_system,L,{Lx,Ly,G{1}},{Cx,Cy,+gradient_amplitudes}, time_grid,rho,'expv-pwc')`.
- Lines 97-98: Select coherence +1; implemented by `rho=coherence(spin_system,rho,{{parameters.spins{1},+1}})`.
- Lines 100-101: Evolve under second positive gradient; implemented by `rho=evolution(spin_system,L+parameters.g_amp*G{1},[],rho,parameters.delta_sml,1,'final')`.
- Lines 103-104: Evolve under delay to obtain complete refocus of the signal; implemented by `rho=step(spin_system,L,rho,parameters.delta)`.
- Lines 106-107: Intensity is the first point in the FID; implemented by `inten=abs(parameters.coil'*rho)`.
- Lines 109-111: Inform the user; implemented by `report(spin_system,['the scaling factor is ' num2str(scaling_factor) ', the gamma B1 max is ' num2str(gamma_B1_max./(2*pi)) ' Hz '])`.

### Key state/data transformations

- Lines 55: computes `L` using `L=H+F+1i*R+1i*K`.
- Lines 58: computes `Lp` using `Lp=operator(spin_system,'L+',parameters.spins{1})`.
- Lines 59: computes `Ly` using `Ly=kron(speye(prod(parameters.npts)),(Lp-Lp')/2i)`.
- Lines 60: computes `Lx` using `Lx=kron(speye(prod(parameters.npts)),(Lp+Lp')/2)`.
- Lines 63: computes `parameters.delta` using `parameters.delta=parameters.delta_big-parameters.delta_sml-parameters.rf_dur`.
- Lines 66: computes `[amps,phases,~,~,scaling_factor]` using `[amps,phases,~,~,scaling_factor]=read_wave(parameters.filename,parameters.pulse_npoints)`.
- Lines 67: computes `time_grid` using `time_grid=parameters.rf_dur*ones(1,parameters.pulse_npoints)/parameters.pulse_npoints`.
- Lines 70: computes `gradient_amplitudes` using `gradient_amplitudes=parameters.sel_g_amp*ones(parameters.pulse_npoints,1)`.
- Lines 74: computes `gamma_B1` using `gamma_B1=parameters.rf_phi/parameters.rf_dur`.
- Lines 75: computes `gamma_B1_max` using `gamma_B1_max=gamma_B1/scaling_factor`.
- Lines 76: computes `amps` using `amps=gamma_B1_max*amps`.
- Lines 79: computes `[Cx,Cy]` using `[Cx,Cy]=polar2cartesian(amps,phases)`.
- Lines 82: computes `rho` using `rho=step(spin_system,Ly,parameters.rho0,pi/2)`.
- Lines 107: computes `inten` using `inten=abs(parameters.coil'*rho)`.

### Local helper functions

- Line 116: `grumble()` — `function grumble(spin_system,parameters,H,R,K,G,F)`.
  - Representative operation: `if ~ismember(spin_system.bas.formalism,{'sphten-liouv'})`.
  - Representative operation: `error('this function is only available for sphten-liouv formalism.')`.

## Outputs

- inten -the absolute value of the first point in
- the free induction decay; this number is
- proportional to the integral of the real
- part of the correctly phased spectrum

## Implementation structure

- A simplified model sequence of the ZS iDOSY pulse sequence. Syntax:
- inten=idosyzs(spin_system,parameters,H,R,K,G,F)
- This sequence must be called from the imaging() context, which
- would provide H, R, K, G and F. Parameters:
- parameters.rho0 -initial state
- parameters.coil -detection state
- parameters.spins -nuclei on which the sequence runs
- parameters.g_amp -gradient amplitude for diffusion
- encoding (T/m)
- parameters.sel_g_amp -gradient amplitude during the
- selective pulse (T/m)
- parameters.rf_phi -phase of the inversion pulse (rad/s)

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `operator()`, `speye()`, `read_wave()`, `polar2cartesian()`, `step()`, `coherence()`, `evolution()`, `shaped_pulse_xy()`, `report()`, `num2str()`, `ismember()`, `ismatrix()`, `all()`, `iscell()`, `isfield()`.
