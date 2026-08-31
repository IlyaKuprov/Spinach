# examples/esr_sol_pulsed/spa_gd_dota_powder.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/esr_sol_pulsed/spa_gd_dota_powder.m`
- Signature: `spa_gd_dota_powder()`
- Total lines: 83

## Purpose

A soft pulse simulation for a gadolinium ion. The soft pulse is simulated using Fokker-Planck formalism. Zero-field split- ting distribution is sampled using the statistical parameters reported in Figure 5 of Raitsimring et al, App. Mag. Res. 28, 281-295 (2005). Powder average simulation with a third-order numerical rotating frame transformation. Calculation time: minutes

## Physical / mathematical content

- Pulsed ESR / EPR solid-state examples. These scripts revolve around electron spin echo sequences, DEER, RIDME, ENDOR, ESEEM, and HYSCORE. They combine anisotropic Zeeman and hyperfine Hamiltonians with selective pulses, echo formation, and orientation averaging.
- The file uses a Fokker-Planck-style enlarged state space in which spatial or orientational coordinates are promoted to extra dimensions and coupled to spin dynamics through differential operators.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 14-15: Preallocate the spectrum; implemented by `spectrum=zeros(2048,1,'like',1i)`.
- Lines 17-18: Get the sampling; implemented by `[D,E,W]=zfs_sampling(30,5,1e-2); drawnow`.
- Lines 20-21: Get the figure going; implemented by `kfigure()`.
- Lines 23-24: Loop over ZFS distribution; implemented by `for n=1:numel(W)`.
- Lines 26-27: Spin system parameters; implemented by `sys.magnet=3.5`.
- Lines 32-33: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 37-38: Disable trajectory-level SSR algorithms; implemented by `sys.disable={'trajlevel'}`.
- Lines 40-41: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 44-45: Sequence parameters; implemented by `parameters.spins={'E8'}`.
- Lines 60-61: Soft pulse parameters; implemented by `parameters.pulse_rnk=2`.
- Lines 68-69: Simulation; implemented by `fid=powder(spin_system,@sp_acquire,parameters,'labframe')`.
- Lines 71-72: Apodisation; implemented by `fid=apodisation(spin_system,fid,{{'exp',10}})`.
- Lines 74-75: Fourier transform and addition; implemented by `spectrum=spectrum+W(n)*fftshift(fft(fid,parameters.zerofill))`.
- Lines 77-78: Plotting; implemented by `plot_1d(spin_system,real(spectrum),parameters); drawnow()`.

### Control flow inferred from the code

- Line 24: `for` loop over `n=1:numel(W)`.

### Key state/data transformations

- Lines 15: computes `spectrum` using `spectrum=zeros(2048,1,'like',1i)`.
- Lines 18: computes `[D,E,W]` using `[D,E,W]=zfs_sampling(30,5,1e-2); drawnow`.
- Lines 27: computes `sys.magnet` using `sys.magnet=3.5`.
- Lines 28: computes `sys.isotopes` using `sys.isotopes={'E8'}`.
- Lines 29: computes `inter.zeeman.scalar` using `inter.zeeman.scalar={2.002319}`.
- Lines 30: computes `inter.coupling.matrix{1,1}` using `inter.coupling.matrix{1,1}=0.56e9*zfs2mat(D(n),E(n),0,0,0)`.
- Lines 33: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 34: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 35: computes `bas.projections` using `bas.projections=-3:3`.
- Lines 38: computes `sys.disable` using `sys.disable={'trajlevel'}`.
- Lines 41: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 45: computes `parameters.spins` using `parameters.spins={'E8'}`.
- Lines 46: computes `parameters.rho0` using `parameters.rho0=state(spin_system,'Lz','E8')`.
- Lines 47: computes `parameters.coil` using `parameters.coil=state(spin_system,'L+','E8')`.
- Lines 48: computes `parameters.decouple` using `parameters.decouple={}`.
- Lines 49: computes `parameters.offset` using `parameters.offset=0`.
- Lines 50: computes `parameters.sweep` using `parameters.sweep=0.8e10`.
- Lines 51: computes `parameters.npoints` using `parameters.npoints=512`.

## Implementation structure

- A soft pulse simulation for a gadolinium ion. The soft pulse
- is simulated using Fokker-Planck formalism. Zero-field split-
- ting distribution is sampled using the statistical parameters
- reported in Figure 5 of Raitsimring et al, App. Mag. Res. 28,
- 281-295 (2005). Powder average simulation with a third-order
- numerical rotating frame transformation.
- Calculation time: minutes
- Preallocate the spectrum
- Get the sampling
- Get the figure going
- Loop over ZFS distribution
- Spin system parameters

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `zfs_sampling()`, `kfigure()`, `zfs2mat()`, `create()`, `basis()`, `state()`, `powder()`, `apodisation()`, `fftshift()`, `plot_1d()`.
