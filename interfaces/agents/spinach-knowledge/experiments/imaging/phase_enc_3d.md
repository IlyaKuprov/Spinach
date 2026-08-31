# experiments/imaging/phase_enc_3d.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/experiments/imaging/phase_enc_3d.m`
- Signature: `fid=phase_enc_3d(spin_system,parameters,H,R,K,G,F)`
- Total lines: 259

## Purpose

3D MRI pulse sequence, with a slice selection stage followed by phase- encoded acquisition of the slice. Syntax: fid=phase_enc_3d(spin_system,parameters,H,R,K,G,F) This sequence must be called from the imaging() context, which would provide H,R,K,G, and F.

## Physical / mathematical content

- Imaging sequence implementations. They build spatially resolved Liouvillians that include gradients, slice-selection RF terms, diffusion, and acquisition operators.
- The file uses a Fokker-Planck-style enlarged state space in which spatial or orientational coordinates are promoted to extra dimensions and coupled to spin dynamics through differential operators.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 41-42: Check consistency; implemented by `grumble(spin_system,parameters,H,R,K,G,F)`.
- Lines 44-45: Assemble the background; implemented by `B=H+F+1i*R+1i*K`.
- Lines 47-48: Make pulse operators and kron them up into Fokker-Planck space; implemented by `Sx=operator(spin_system,'Lx','1H')`.
- Lines 51-52: Make and inflate polyadics; implemented by `spc_dim=prod(parameters.npts)`.
- Lines 59-60: Slice selection operator; implemented by `L=B+parameters.ss_grad_amp*G{1}`.
- Lines 62-65: Apply the slice selection pulse; implemented by `parameters.rho0=shaped_pulse_af(spin_system,L,Sx,Sy,parameters.rho0, parameters.rf_frq_list,parameters.rf_amp_list, parameters.rf_dur_list,parameters.rf_phi,2,'expv')`.
- Lines 67-68: Rollback operator; implemented by `L=B-parameters.ss_grad_amp*G{1}`.
- Lines 70-72: Run the rollback gradient; implemented by `parameters.rho0=evolution(spin_system,L,[],parameters.rho0, sum(parameters.rf_dur_list)/2,1,'final')`.
- Lines 74-75: Evolve for the echo time; implemented by `parameters.rho0=evolution(spin_system,B,[],parameters.rho0,parameters.t_echo,1,'final')`.
- Lines 77-78: Apply 180-degree pulse; implemented by `parameters.rho0=step(spin_system,Sy,parameters.rho0,pi)`.
- Lines 83-84: Project out Hx+1i*Hy in every voxel; implemented by `mri_slice=fpl2phan(parameters.rho0,state(spin_system,'L+','1H'),parameters.npts)`.
- Lines 86-87: Get sample dimension information; implemented by `dims=zeros(1,6)`.
- Lines 91-92: Draw the slice in three dimensions; implemented by `kfigure(); volplot(abs(mri_slice),dims)`.
- Lines 95-98: Get phase encoding gradient range; implemented by `pe_grad_amps=linspace(-parameters.pe_grad_amp, parameters.pe_grad_amp, parameters.image_size(1))`.
- Lines 100-101: Get k-space sampling parameters; implemented by `nsteps=parameters.image_size(2)-1; step_length=parameters.ro_grad_dur/nsteps`.
- Lines 103-104: Preallocate the image; implemented by `fid=zeros(parameters.image_size,'like',1i)`.
- Lines 106-107: Precompute evolution generators; implemented by `F_preroll=B-parameters.ro_grad_amp*G{3}`.
- Lines 110-111: Loop over frequencies; implemented by `parfor n=1:parameters.image_size(1)`.

### Control flow inferred from the code

- Line 55: conditional branch on `~ismember('polyadic',spin_system.sys.enable)`.
- Line 111: `parfor` loop over `n=1:parameters.image_size(1)`.

### Key state/data transformations

- Lines 45: computes `B` using `B=H+F+1i*R+1i*K`.
- Lines 48: computes `Sx` using `Sx=operator(spin_system,'Lx','1H')`.
- Lines 49: computes `Sy` using `Sy=operator(spin_system,'Ly','1H')`.
- Lines 52: computes `spc_dim` using `spc_dim=prod(parameters.npts)`.
- Lines 60: computes `L` using `L=B+parameters.ss_grad_amp*G{1}`.
- Lines 63-65: computes `parameters.rho0` using `parameters.rho0=shaped_pulse_af(spin_system,L,Sx,Sy,parameters.rho0, parameters.rf_frq_list,parameters.rf_amp_list, parameters.rf_dur_list,parameters.rf_phi,2,'expv')`.
- Lines 84: computes `mri_slice` using `mri_slice=fpl2phan(parameters.rho0,state(spin_system,'L+','1H'),parameters.npts)`.
- Lines 87: computes `dims` using `dims=zeros(1,6)`.
- Lines 88: computes `dims([1 3 5])` using `dims([1 3 5])=-parameters.dims/2`.
- Lines 89: computes `dims([2 4 6])` using `dims([2 4 6])=+parameters.dims/2`.
- Lines 96-98: computes `pe_grad_amps` using `pe_grad_amps=linspace(-parameters.pe_grad_amp, parameters.pe_grad_amp, parameters.image_size(1))`.
- Lines 101: computes `nsteps` using `nsteps=parameters.image_size(2)-1; step_length=parameters.ro_grad_dur/nsteps`.
- Lines 104: computes `fid` using `fid=zeros(parameters.image_size,'like',1i)`.
- Lines 107: computes `F_preroll` using `F_preroll=B-parameters.ro_grad_amp*G{3}`.
- Lines 108: computes `F_readout` using `F_readout=B+parameters.ro_grad_amp*G{3}`.
- Lines 114: computes `rho` using `rho=evolution(spin_system,B+pe_grad_amps(n)*G{2},[],parameters.rho0,parameters.pe_grad_dur,1,'final')`.
- Lines 120: computes `fid(n,:)` using `fid(n,:)=evolution(spin_system,F_readout,parameters.coil,rho,step_length,nsteps,'observable')`.

### Local helper functions

- Line 127: `grumble()` — `function grumble(spin_system,parameters,H,R,K,G,F)`.
  - Representative operation: `if ~ismember(spin_system.bas.formalism,{'sphten-liouv','zeeman-liouv'})`.
  - Representative operation: `error('this function is only available in sphten-liouv and zeeman-liouv formalisms.')`.

## Parameters / inputs

- parameters.ss_grad_amp -the amplitude of slice selection
- gradient,T/m
- parameters.pe_grad_amp -phase encoding gradient amplitude, T/m
- parameters.ro_grad_amp -readout gradient amplitude, T/m
- parameters.ss_grad_dur -the duration of the slice selection
- gradient, seconds
- parameters.pe_grad_dur -the duration of the phase encoding
- gradient, seconds
- parameters.ro_grad_dur -the duration of the readout gradient,
- seconds
- parameters.image_size -number of points in each dimension of
- the resulting image

## Outputs

- fid -k-space representation of the image

## Implementation structure

- 3D MRI pulse sequence, with a slice selection stage followed by phase-
- encoded acquisition of the slice. Syntax:
- fid=phase_enc_3d(spin_system,parameters,H,R,K,G,F)
- This sequence must be called from the imaging() context, which would
- provide H,R,K,G, and F.
- parameters.ss_grad_amp - the amplitude of slice selection
- gradient,T/m
- parameters.pe_grad_amp - phase encoding gradient amplitude, T/m
- parameters.ro_grad_amp - readout gradient amplitude, T/m
- parameters.ss_grad_dur - the duration of the slice selection
- gradient, seconds
- parameters.pe_grad_dur - the duration of the phase encoding

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `operator()`, `polyadic()`, `opium()`, `ismember()`, `inflate()`, `shaped_pulse_af()`, `evolution()`, `step()`, `fpl2phan()`, `state()`, `dims()`, `kfigure()`, `volplot()`, `ktitle()`, `pe_grad_amps()`.
