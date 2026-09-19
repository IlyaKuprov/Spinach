# experiments/imaging/epi_3d.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/experiments/imaging/epi_3d.m`
- Signature: `fid=epi_3d(spin_system,parameters,H,R,K,G,F)`
- Total lines: 348

## Purpose

Diffusion weighted 3D echo planar imaging pulse sequence. Syntax: fid=epi_3d(spin_system,parameters,H,R,K,G,F) This sequence must be called from the imaging() context, which would provide H,R,K,G, and F. Parameters: parameters.image_size -number of points in each dimension of the resulting image parameters.ss_grad_amp -the amplitude of slice selection gradient,T/m parameters.pe_grad_amp -phase encoding gradient ampli

## Physical / mathematical content

- Imaging sequence implementations. They build spatially resolved Liouvillians that include gradients, slice-selection RF terms, diffusion, and acquisition operators.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 40-41: Check consistency; implemented by `grumble(spin_system,parameters,H,R,K,G,F)`.
- Lines 43-44: Assemble the background; implemented by `B=H+F+1i*R+1i*K`.
- Lines 46-47: Make pulse operators; implemented by `Hx=operator(spin_system,'Lx','1H')`.
- Lines 52-56: Apply the slice selection pulse; implemented by `rho=shaped_pulse_af(spin_system,B+parameters.ss_grad_amp*G{1}, Hx,Hy,parameters.rho0,parameters.rf_frq_list, parameters.rf_amp_list,parameters.rf_dur_list, parameters.rf…`.
- Lines 58-60: Run the rollback gradient; implemented by `rho=evolution(spin_system,B-parameters.ss_grad_amp*G{1}, [],rho,sum(parameters.rf_dur_list)/2,1,'final')`.
- Lines 62-63: Optional diffusion gradient; implemented by `if isfield(parameters,'diff_g_amp')`.
- Lines 65-69: Evolve for the echo time in the presence of diffusion gradient; implemented by `rho=step(spin_system,B+parameters.diff_g_amp(1)*G{1}+ parameters.diff_g_amp(2)*G{2}+ parameters.diff_g_amp(3)*G{3}, rho,parameters.t_echo)`.
- Lines 73-74: Simply evolve for the echo time; implemented by `rho=evolution(spin_system,B,[],rho,parameters.t_echo,1,'final')`.
- Lines 78-79: Apply an ideal 180-degree pulse; implemented by `rho=step(spin_system,Hy,rho,pi)`.
- Lines 97-98: Project out Hx+1i*Hy in every voxel; implemented by `mri_slice=fpl2phan(rho,state(spin_system,'L+','1H'),parameters.npts)`.
- Lines 100-101: Get sample dimension information; implemented by `dims=zeros(1,6)`.
- Lines 105-106: Draw the slice in three dimensions; implemented by `kfigure(); volplot(abs(mri_slice),dims)`.
- Lines 109-112: Preroll both gradients together for the duration of the shorter prephaser; implemented by `rho=evolution(spin_system,B-parameters.pe_grad_amp*G{2} -parameters.ro_grad_amp*G{3},[],rho, min(parameters.pe_grad_dur,parameters.ro_grad_dur)/2,1,'final')`.
- Lines 114-121: Finish the longer prephaser on its own; implemented by `if parameters.pe_grad_dur>parameters.ro_grad_dur`.
- Lines 123-125: Precompute propagators; implemented by `P_ro_p=propagator(spin_system,B+parameters.ro_grad_amp*G{3}, parameters.ro_grad_dur/(parameters.image_size(2)-1))`.
- Lines 131-132: Move to GPU if necessary; implemented by `if ismember('gpu',spin_system.sys.enable)`.
- Lines 134-135: Upload relevant objects; implemented by `rho=gpuArray(rho); P_ro_p=gpuArray(P_ro_p)`.
- Lines 139-140: Preallocate k-space image; implemented by `fid=1i*gpuArray.zeros(parameters.image_size)`.
- Lines 142-143: Inform the user; implemented by `report(spin_system,'running the EPI loop on GPU ')`.

### Control flow inferred from the code

- Line 115: conditional branch on `parameters.pe_grad_dur>parameters.ro_grad_dur`, with an `elseif` for `parameters.ro_grad_dur>parameters.pe_grad_dur`.
- Line 63: conditional branch on `isfield(parameters,'diff_g_amp')`.
- Line 63: conditional branch on `isfield(parameters,'diff_g_amp')`.
- Line 132: conditional branch on `ismember('gpu',spin_system.sys.enable)`.
- Line 159: `for` loop over `n=1:parameters.image_size(1)`.
- Line 169: `for` loop over `k=1:parameters.image_size(2)`.
- Line 172: conditional branch on `ro_grad_sign>0`.

### Key state/data transformations

- Lines 44: computes `B` using `B=H+F+1i*R+1i*K`.
- Lines 47: computes `Hx` using `Hx=operator(spin_system,'Lx','1H')`.
- Lines 48: computes `Hy` using `Hy=operator(spin_system,'Ly','1H')`.
- Lines 53-56: computes `rho` using `rho=shaped_pulse_af(spin_system,B+parameters.ss_grad_amp*G{1}, Hx,Hy,parameters.rho0,parameters.rf_frq_list, parameters.rf_amp_list,parameters.rf_dur_list, parameters.rf…`.
- Lines 98: computes `mri_slice` using `mri_slice=fpl2phan(rho,state(spin_system,'L+','1H'),parameters.npts)`.
- Lines 101: computes `dims` using `dims=zeros(1,6)`.
- Lines 102: computes `dims([1 3 5])` using `dims([1 3 5])=-parameters.dims/2`.
- Lines 103: computes `dims([2 4 6])` using `dims([2 4 6])=+parameters.dims/2`.
- Lines 124-125: computes `P_ro_p` using `P_ro_p=propagator(spin_system,B+parameters.ro_grad_amp*G{3}, parameters.ro_grad_dur/(parameters.image_size(2)-1))`.
- Lines 126-127: computes `P_ro_m` using `P_ro_m=propagator(spin_system,B-parameters.ro_grad_amp*G{3}, parameters.ro_grad_dur/(parameters.image_size(2)-1))`.
- Lines 128-129: computes `P_pe_p` using `P_pe_p=propagator(spin_system,B+parameters.pe_grad_amp*G{2}, parameters.pe_grad_dur/(parameters.image_size(1)-1))`.
- Lines 137: computes `coil` using `coil=gpuArray(parameters.coil)`.
- Lines 140: computes `fid` using `fid=1i*gpuArray.zeros(parameters.image_size)`.
- Lines 166: computes `ro_grad_sign` using `ro_grad_sign=2*mod(n,2)-1`.
- Lines 173: computes `fid(n,k)` using `fid(n,k)=hdot(coil,rho)`.
- Lines 177: computes `fid(n,parameters.image_size(2)-k+1)` using `fid(n,parameters.image_size(2)-k+1)=hdot(coil,rho)`.

### Local helper functions

- Line 193: `grumble()` — `function grumble(spin_system,parameters,H,R,K,G,F)`.
  - Representative operation: `if ~ismember(spin_system.bas.formalism,{'sphten-liouv','zeeman-liouv'})`.
  - Representative operation: `error('this function is only available in sphten-liouv and zeeman-liouv formalisms.')`.

## Outputs

- fid -k-space representation of the image

## Implementation structure

- Diffusion weighted 3D echo planar imaging pulse sequence. Syntax:
- fid=epi_3d(spin_system,parameters,H,R,K,G,F)
- This sequence must be called from the imaging() context, which
- would provide H,R,K,G, and F. Parameters:
- parameters.image_size - number of points in each dimension of
- the resulting image
- parameters.ss_grad_amp - the amplitude of slice selection
- gradient,T/m
- parameters.pe_grad_amp - phase encoding gradient amplitude, T/m
- parameters.pe_grad_dur - phase encoding gradient duration, s
- parameters.ro_grad_amp - readout gradient amplitude, T/m
- parameters.ro_grad_dur - readout gradient duration, s

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `operator()`, `speye()`, `shaped_pulse_af()`, `evolution()`, `isfield()`, `step()`, `fpl2phan()`, `state()`, `dims()`, `kfigure()`, `volplot()`, `ktitle()`, `propagator()`, `ismember()`, `gpuArray()`.
