# examples/esr_sol_pulsed/hpa_triplet.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/esr_sol_pulsed/hpa_triplet.m`
- Signature: `hpa_triplet()`
- Total lines: 73

## Purpose

Hypothetical powder averaged X-band pulse-acquire ESR spectrum of photogenerated pentacene triplet state. Calculation time: seconds.

## Physical / mathematical content

- Pulsed ESR / EPR solid-state examples. These scripts revolve around electron spin echo sequences, DEER, RIDME, ENDOR, ESEEM, and HYSCORE. They combine anisotropic Zeeman and hyperfine Hamiltonians with selective pulses, echo formation, and orientation averaging.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 11-12: Magnet field; implemented by `sys.magnet=0.33`.
- Lines 14-15: Triplet electron; implemented by `sys.isotopes={'E3'}`.
- Lines 17-18: Zeeman tensor, assumed isotropic; implemented by `inter.zeeman.matrix={diag([2.0 2.0 2.0])}`.
- Lines 20-21: ZFS, photo-excited pentacene triplet; implemented by `D=1360.1*1e6; E=-47.2*1e6`.
- Lines 24-25: Basis set; implemented by `bas.formalism='zeeman-hilb'`.
- Lines 28-29: Disable trajectory-level SSR algorithms; implemented by `sys.disable={'trajlevel'}`.
- Lines 31-32: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 35-36: Sequence parameters; implemented by `parameters.spins={'E3'}`.
- Lines 49-50: Zeeman tensor into Hz/Tesla; implemented by `Z=-spin('E')*inter.zeeman.matrix{1,1}/(2*pi*2.00231930436256)`.
- Lines 52-58: Orientation-dependent initial condition; implemented by `parameters.rho0=@(alp,bet,gam)zftrip(spin_system,euler2dcm(alp,bet,gam)* inter.coupling.matrix{1,1}* euler2dcm(alp,bet,gam)', [0.56 0.31 0.13], euler2dcm(alp,bet,gam)*Z*…`.
- Lines 60-61: Simulation of a pulse-acquire experiment; implemented by `fid=powder(spin_system,@hp_acquire,parameters,'esr')`.
- Lines 63-64: Apodisation; implemented by `fid=apodisation(spin_system,fid,{{'crisp'}})`.
- Lines 66-67: Fourier transform; implemented by `spectrum=fftshift(fft(fid,parameters.zerofill))`.
- Lines 69-70: Plotting; implemented by `kfigure(); plot_1d(spin_system,real(spectrum),parameters)`.

### Key state/data transformations

- Lines 12: computes `sys.magnet` using `sys.magnet=0.33`.
- Lines 15: computes `sys.isotopes` using `sys.isotopes={'E3'}`.
- Lines 18: computes `inter.zeeman.matrix` using `inter.zeeman.matrix={diag([2.0 2.0 2.0])}`.
- Lines 21: computes `D` using `D=1360.1*1e6; E=-47.2*1e6`.
- Lines 22: computes `inter.coupling.matrix` using `inter.coupling.matrix={zfs2mat(D,E,0,0,0)}`.
- Lines 25: computes `bas.formalism` using `bas.formalism='zeeman-hilb'`.
- Lines 26: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 29: computes `sys.disable` using `sys.disable={'trajlevel'}`.
- Lines 32: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 36: computes `parameters.spins` using `parameters.spins={'E3'}`.
- Lines 37: computes `parameters.coil` using `parameters.coil=state(spin_system,'L+','E3')`.
- Lines 38: computes `parameters.pulse_op` using `parameters.pulse_op=operator(spin_system,'Ly','E3')`.
- Lines 39: computes `parameters.pulse_angle` using `parameters.pulse_angle=pi/4`.
- Lines 40: computes `parameters.offset` using `parameters.offset=0`.
- Lines 41: computes `parameters.sweep` using `parameters.sweep=4e9`.
- Lines 42: computes `parameters.npoints` using `parameters.npoints=128`.
- Lines 43: computes `parameters.zerofill` using `parameters.zerofill=512`.
- Lines 44: computes `parameters.axis_units` using `parameters.axis_units='GHz-labframe'`.

## Implementation structure

- Hypothetical powder averaged X-band pulse-acquire ESR
- spectrum of photogenerated pentacene triplet state.
- Calculation time: seconds.
- Magnet field
- Triplet electron
- Zeeman tensor, assumed isotropic
- ZFS, photo-excited pentacene triplet
- Basis set
- Disable trajectory-level SSR algorithms
- Spinach housekeeping
- Sequence parameters
- Zeeman tensor into Hz/Tesla

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `zfs2mat()`, `create()`, `basis()`, `state()`, `operator()`, `spin()`, `zftrip()`, `euler2dcm()`, `powder()`, `apodisation()`, `fftshift()`, `kfigure()`, `plot_1d()`.
