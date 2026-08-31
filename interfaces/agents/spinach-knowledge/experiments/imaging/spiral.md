# experiments/imaging/spiral.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/experiments/imaging/spiral.m`
- Signature: `mri=spiral(spin_system,parameters,H,R,K,G,F)`
- Total lines: 211

## Purpose

2D imaging sequence with spiral sampling of the k-space. Syntax: mri=spiral_pulse_sequence(spin_system,parameters,H,R,K,G,F) This sequence must be called from the imaging() context, which would provide H,R,K,G, and F.

## Physical / mathematical content

- Imaging sequence implementations. They build spatially resolved Liouvillians that include gradients, slice-selection RF terms, diffusion, and acquisition operators.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 35-36: Check consistency; implemented by `grumble(spin_system,parameters,H,R,K,G,F)`.
- Lines 38-39: Compose Liouvillian; implemented by `L=H+F+1i*R+1i*K`.
- Lines 41-42: Make pulse operators; implemented by `Hp=operator(spin_system,'L+',parameters.spins{1})`.
- Lines 45-46: Apply 90-degree pulse; implemented by `rho=step(spin_system,Hy,parameters.rho0,pi/2)`.
- Lines 48-49: Evolve the for the echo time; implemented by `rho=evolution(spin_system,L,[],rho,parameters.t_echo,1,'final')`.
- Lines 51-52: Apply 180-degree pulse; implemented by `rho=step(spin_system,Hy,rho,pi)`.
- Lines 57-58: Build the spiral; implemented by `time_grid=(0:(parameters.spiral_npts-1))*(parameters.spiral_dur/parameters.spiral_npts)`.
- Lines 62-63: Get the time step; implemented by `time_step=(parameters.spiral_dur/parameters.spiral_npts)`.
- Lines 65-66: Accumulate the k-space trajectory in cycles per metre; implemented by `gam=spin(parameters.spins{1})`.
- Lines 70-71: Report the sampling spiral to the user; implemented by `kfigure(); plot(spiral_kx,spiral_ky,'b-'); ktitle('$k$-space sampling spiral')`.
- Lines 75-76: Prelocate the spiral trajectory array; implemented by `spiral_z=zeros([parameters.spiral_npts 1],'like',1i)`.
- Lines 78-79: Start feedback timer; implemented by `feedback=tic()`.
- Lines 81-82: Propagate through the spiral; implemented by `for n=1:parameters.spiral_npts`.
- Lines 84-85: Record trajectory point; implemented by `spiral_z(n)=parameters.coil'*rho`.
- Lines 87-88: Take a step forward; implemented by `rho=step(spin_system,L+spiral_x(n)*G{1}+spiral_y(n)*G{2},rho,time_step)`.
- Lines 90-91: Inform the user; implemented by `if (n==parameters.spiral_npts)||(toc(feedback)>1)`.
- Lines 99-100: Resample onto square grid; implemented by `k_max=max(hypot(spiral_kx,spiral_ky))`.
- Lines 105-106: Replace NaN values with zeros; implemented by `fid(isnan(fid))=0`.

### Control flow inferred from the code

- Line 82: `for` loop over `n=1:parameters.spiral_npts`.
- Line 91: conditional branch on `(n==parameters.spiral_npts)||(toc(feedback)>1)`.

### Key state/data transformations

- Lines 39: computes `L` using `L=H+F+1i*R+1i*K`.
- Lines 42: computes `Hp` using `Hp=operator(spin_system,'L+',parameters.spins{1})`.
- Lines 43: computes `Hy` using `Hy=kron(speye(prod(parameters.npts)),(Hp-Hp')/2i)`.
- Lines 46: computes `rho` using `rho=step(spin_system,Hy,parameters.rho0,pi/2)`.
- Lines 58: computes `time_grid` using `time_grid=(0:(parameters.spiral_npts-1))*(parameters.spiral_dur/parameters.spiral_npts)`.
- Lines 59: computes `spiral_x` using `spiral_x=time_grid*(parameters.grad_amp/parameters.spiral_dur).*cos(parameters.spiral_frq*time_grid)`.
- Lines 60: computes `spiral_y` using `spiral_y=time_grid*(parameters.grad_amp/parameters.spiral_dur).*sin(parameters.spiral_frq*time_grid)`.
- Lines 63: computes `time_step` using `time_step=(parameters.spiral_dur/parameters.spiral_npts)`.
- Lines 66: computes `gam` using `gam=spin(parameters.spins{1})`.
- Lines 67: computes `spiral_kx` using `spiral_kx=(gam/(2*pi))*time_step*[0 cumsum(spiral_x(1:(end-1)))]`.
- Lines 68: computes `spiral_ky` using `spiral_ky=(gam/(2*pi))*time_step*[0 cumsum(spiral_y(1:(end-1)))]`.
- Lines 76: computes `spiral_z` using `spiral_z=zeros([parameters.spiral_npts 1],'like',1i)`.
- Lines 79: computes `feedback` using `feedback=tic()`.
- Lines 85: computes `spiral_z(n)` using `spiral_z(n)=parameters.coil'*rho`.
- Lines 100: computes `k_max` using `k_max=max(hypot(spiral_kx,spiral_ky))`.
- Lines 101: computes `x_grid` using `x_grid=linspace(-k_max,k_max,sqrt(parameters.spiral_npts)/2)`.
- Lines 102: computes `y_grid` using `y_grid=linspace(-k_max,k_max,sqrt(parameters.spiral_npts)/2)`.
- Lines 103: computes `[X,Y]` using `[X,Y]=ndgrid(x_grid,y_grid); fid=griddata(spiral_kx,spiral_ky,spiral_z,X,Y,'cubic')`.

### Local helper functions

- Line 117: `grumble()` — `function grumble(spin_system,parameters,H,R,K,G,F)`.
  - Representative operation: `if ~ismember(spin_system.bas.formalism,{'sphten-liouv','zeeman-liouv'})`.
  - Representative operation: `error('this function is only available in sphten-liouv and zeeman-liouv formalisms.')`.

## Parameters / inputs

- parameters.t_echo -echo time, seconds
- parameters.spiral_frq -frequency of the spiral, Hz
- parameters.spiral_dur -duration of the spiral, seconds
- parameters.grad_amp -gradient amplitude at the end
- of the spiral, T/m

## Outputs

- mri -MRI image with square sinebell apodisation
- Prior to being Fourier transformed, the spiral is resampled onto
- a suitable square k-space grid -change the sequence to output
- spiral_kx, spiral_ky and spiral_z variables if you plan to pro-
- cess the data in a different way.

## Implementation structure

- 2D imaging sequence with spiral sampling of the k-space. Syntax:
- mri=spiral_pulse_sequence(spin_system,parameters,H,R,K,G,F)
- This sequence must be called from the imaging() context, which
- would provide H,R,K,G, and F.
- parameters.t_echo -echo time, seconds
- parameters.spiral_frq -frequency of the spiral, Hz
- parameters.spiral_dur -duration of the spiral, seconds
- parameters.grad_amp -gradient amplitude at the end
- of the spiral, T/m
- mri -MRI image with square sinebell apodisation
- Prior to being Fourier transformed, the spiral is resampled onto
- a suitable square k-space grid -change the sequence to output

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `operator()`, `speye()`, `step()`, `evolution()`, `spin()`, `cumsum()`, `spiral_x()`, `spiral_y()`, `kfigure()`, `ktitle()`, `kxlabel()`, `kylabel()`, `tic()`, `spiral_z()`, `toc()`.
