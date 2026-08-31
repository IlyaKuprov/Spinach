# examples/esr_sol_pulsed/holeburn_gd_dota_powder.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/esr_sol_pulsed/holeburn_gd_dota_powder.m`
- Signature: `holeburn_gd_dota_powder()`
- Total lines: 93

## Purpose

A hole burning simulation for a gadolinium ion. The soft pulse is simulated using Fokker-Planck formalism. Zero-field splitting dis- ribution is sampled using the statistical parameters reported in Figure 5 of Raitsimring et al, App. Mag. Res. 28, 281-295 (2005). A numerical powder grid and numerical second-order rotating frame transformation are used. Note: non-central transition Gd(III) holes are very shallow. Calc

## Physical / mathematical content

- Pulsed ESR / EPR solid-state examples. These scripts revolve around electron spin echo sequences, DEER, RIDME, ENDOR, ESEEM, and HYSCORE. They combine anisotropic Zeeman and hyperfine Hamiltonians with selective pulses, echo formation, and orientation averaging.
- The file uses a Fokker-Planck-style enlarged state space in which spatial or orientational coordinates are promoted to extra dimensions and coupled to spin dynamics through differential operators.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 16-17: Initialize the spectra; implemented by `spectrum_a=zeros(2048,1,'like',1i)`.
- Lines 20-21: Get the sampling; implemented by `[D,E,W]=zfs_sampling(30,5,1e-2); drawnow()`.
- Lines 23-24: Get the figure going; implemented by `kfigure()`.
- Lines 26-27: Loop over ZFS distribution; implemented by `for n=1:numel(W)`.
- Lines 29-30: Spin system parameters; implemented by `sys.magnet=3.5`.
- Lines 35-36: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 40-41: Disable trajectory-level SSR algorithms; implemented by `sys.disable={'trajlevel'}`.
- Lines 43-44: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 47-48: Sequence parameters; implemented by `parameters.spins={'E8'}`.
- Lines 63-64: Soft pulse parameters; implemented by `parameters.pulse_rnk=2`.
- Lines 70-71: Simulation A; implemented by `parameters.pulse_pwr=0`.
- Lines 74-75: Simulation B; implemented by `parameters.pulse_pwr=2*pi*1e7`.
- Lines 78-79: Apodisation; implemented by `fid_a=apodisation(spin_system,fid_a,{{'exp',10}})`.
- Lines 82-83: Fourier transform and addition; implemented by `spectrum_a=spectrum_a+W(n)*fftshift(fft(fid_a,parameters.zerofill))`.
- Lines 86-87: Plotting; implemented by `plot_1d(spin_system,real(spectrum_a),parameters,'r-'); hold on`.

### Control flow inferred from the code

- Line 27: `for` loop over `n=1:numel(W)`.

### Key state/data transformations

- Lines 17: computes `spectrum_a` using `spectrum_a=zeros(2048,1,'like',1i)`.
- Lines 18: computes `spectrum_b` using `spectrum_b=zeros(2048,1,'like',1i)`.
- Lines 21: computes `[D,E,W]` using `[D,E,W]=zfs_sampling(30,5,1e-2); drawnow()`.
- Lines 30: computes `sys.magnet` using `sys.magnet=3.5`.
- Lines 31: computes `sys.isotopes` using `sys.isotopes={'E8'}`.
- Lines 32: computes `inter.zeeman.scalar` using `inter.zeeman.scalar={2.002319}`.
- Lines 33: computes `inter.coupling.matrix{1,1}` using `inter.coupling.matrix{1,1}=0.56e9*zfs2mat(D(n),E(n),0,0,0)`.
- Lines 36: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 37: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 38: computes `bas.projections` using `bas.projections=-3:3`.
- Lines 41: computes `sys.disable` using `sys.disable={'trajlevel'}`.
- Lines 44: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 48: computes `parameters.spins` using `parameters.spins={'E8'}`.
- Lines 49: computes `parameters.rho0` using `parameters.rho0=state(spin_system,'Lz','E8')`.
- Lines 50: computes `parameters.coil` using `parameters.coil=state(spin_system,'L+','E8')`.
- Lines 51: computes `parameters.decouple` using `parameters.decouple={}`.
- Lines 52: computes `parameters.offset` using `parameters.offset=0`.
- Lines 53: computes `parameters.sweep` using `parameters.sweep=0.8e10`.

## Implementation structure

- A hole burning simulation for a gadolinium ion. The soft pulse is
- simulated using Fokker-Planck formalism. Zero-field splitting dis-
- ribution is sampled using the statistical parameters reported in
- Figure 5 of Raitsimring et al, App. Mag. Res. 28, 281-295 (2005).
- A numerical powder grid and numerical second-order rotating frame
- transformation are used.
- Note: non-central transition Gd(III) holes are very shallow.
- Calculation time: minutes
- Initialize the spectra
- Get the sampling
- Get the figure going
- Loop over ZFS distribution

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `zfs_sampling()`, `kfigure()`, `zfs2mat()`, `create()`, `basis()`, `state()`, `powder()`, `apodisation()`, `fftshift()`, `plot_1d()`.
