# examples/shaped_pulses/shaped_pulse_chirp_xy.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/shaped_pulses/shaped_pulse_chirp_xy.m`
- Signature: `shaped_pulse_chirp_xy()`
- Total lines: 86

## Purpose

Chirped inversion pulse. Calculation time: seconds

## Physical / mathematical content

- Shaped-pulse examples. These scripts demonstrate amplitude, phase, frequency, and gradient waveform design, including adiabatic sweeps, excitation profiles, and hardware-response considerations.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 9-10: Magnetic field; implemented by `sys.magnet=14.1`.
- Lines 12-15: Isotopes; implemented by `sys.isotopes={'1H','1H','1H','1H','1H','1H','1H','1H','1H','1H','1H', '1H','1H','1H','1H','1H','1H','1H','1H','1H','1H','1H', '1H','1H','1H','1H','1H','1H','1H','1H','1H…`.
- Lines 17-18: Zeeman interactions; implemented by `inter.zeeman.scalar=num2cell(linspace(-4,4,31))`.
- Lines 20-21: Couplings; implemented by `inter.coupling.scalar=cell(31)`.
- Lines 26-27: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 32-33: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 36-37: Set up acquisition; implemented by `parameters.spins={'1H'}`.
- Lines 47-48: Pulse infrastructure; implemented by `H=hamiltonian(assume(spin_system,'nmr'))`.
- Lines 54-55: Chirp waveform in amplitude-frequency coordinates; implemented by `[Cx,Cy,durs]=chirp_pulse(500,0.1,2000,20,'wurst-adaptive')`.
- Lines 57-59: Soft pulse; implemented by `parameters.rho0=shaped_pulse_xy(spin_system,H,{Lx,Ly},{Cx,Cy}, durs,parameters.rho0,'expv-pwc')`.
- Lines 61-62: Homospoil; implemented by `parameters.rho0=homospoil(spin_system,parameters.rho0,'destroy')`.
- Lines 64-65: Global hard pulse; implemented by `parameters.rho0=step(spin_system,Ly,parameters.rho0,pi/2)`.
- Lines 67-68: Acquisition; implemented by `fid=acquire(spin_system,parameters,H,R,K)`.
- Lines 70-71: Apodisation; implemented by `fid=apodisation(spin_system,fid,{{'exp',6}})`.
- Lines 73-74: Fourier transform; implemented by `spectrum=fftshift(fft(fid,parameters.zerofill))`.
- Lines 77-78: Plotting; implemented by `kfigure(); scale_figure([2.0 1.0])`.

### Control flow inferred from the code

- Line 22: `for` loop over `n=1:30`.

### Key state/data transformations

- Lines 10: computes `sys.magnet` using `sys.magnet=14.1`.
- Lines 13-15: computes `sys.isotopes` using `sys.isotopes={'1H','1H','1H','1H','1H','1H','1H','1H','1H','1H','1H', '1H','1H','1H','1H','1H','1H','1H','1H','1H','1H','1H', '1H','1H','1H','1H','1H','1H','1H','1H','1H…`.
- Lines 18: computes `inter.zeeman.scalar` using `inter.zeeman.scalar=num2cell(linspace(-4,4,31))`.
- Lines 21: computes `inter.coupling.scalar` using `inter.coupling.scalar=cell(31)`.
- Lines 23: computes `inter.coupling.scalar{n,n+1}` using `inter.coupling.scalar{n,n+1}=20`.
- Lines 27: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 28: computes `bas.approximation` using `bas.approximation='IK-2'`.
- Lines 29: computes `bas.connectivity` using `bas.connectivity='scalar_couplings'`.
- Lines 30: computes `bas.space_level` using `bas.space_level=1`.
- Lines 33: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 37: computes `parameters.spins` using `parameters.spins={'1H'}`.
- Lines 38: computes `parameters.rho0` using `parameters.rho0=state(spin_system,'Lz','1H')`.
- Lines 39: computes `parameters.coil` using `parameters.coil=state(spin_system,'L+','1H')`.
- Lines 40: computes `parameters.decouple` using `parameters.decouple={}`.
- Lines 41: computes `parameters.offset` using `parameters.offset=0`.
- Lines 42: computes `parameters.sweep` using `parameters.sweep=5100`.
- Lines 43: computes `parameters.npoints` using `parameters.npoints=2048`.
- Lines 44: computes `parameters.zerofill` using `parameters.zerofill=8192`.

## Implementation structure

- Chirped inversion pulse.
- Calculation time: seconds
- Magnetic field
- Isotopes
- Zeeman interactions
- Couplings
- Basis set
- Spinach housekeeping
- Set up acquisition
- Pulse infrastructure
- Chirp waveform in amplitude-frequency coordinates
- Soft pulse

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `num2cell()`, `create()`, `basis()`, `state()`, `hamiltonian()`, `assume()`, `relaxation()`, `kinetics()`, `operator()`, `chirp_pulse()`, `shaped_pulse_xy()`, `homospoil()`, `step()`, `acquire()`, `apodisation()`, `fftshift()`.
