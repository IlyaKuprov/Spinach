# examples/shaped_pulses/shaped_pulse_chirp_af.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/shaped_pulses/shaped_pulse_chirp_af.m`
- Signature: `shaped_pulse_chirp_af()`
- Total lines: 87

## Purpose

Chirp inversion pulse using the Fokker-Planck formalism. Fewer points are required by the amplitude-frequency method than the "two points per period of the largest frequency" Nyquist-Shan- non condition would need for the {Cx,Cy} parameterised simula- tion of a chirped pulse. Calculation time: seconds

## Physical / mathematical content

- Shaped-pulse examples. These scripts demonstrate amplitude, phase, frequency, and gradient waveform design, including adiabatic sweeps, excitation profiles, and hardware-response considerations.
- The file uses a Fokker-Planck-style enlarged state space in which spatial or orientational coordinates are promoted to extra dimensions and coupled to spin dynamics through differential operators.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 13-14: Magnetic field; implemented by `sys.magnet=14.1`.
- Lines 16-19: Isotopes; implemented by `sys.isotopes={'1H','1H','1H','1H','1H','1H','1H','1H','1H','1H','1H', '1H','1H','1H','1H','1H','1H','1H','1H','1H','1H','1H', '1H','1H','1H','1H','1H','1H','1H','1H','1H…`.
- Lines 21-22: Zeeman interactions; implemented by `inter.zeeman.scalar=num2cell(linspace(-4,4,31))`.
- Lines 24-25: Couplings; implemented by `inter.coupling.scalar=cell(31)`.
- Lines 30-31: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 36-37: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 40-41: Set up acquisition; implemented by `parameters.spins={'1H'}`.
- Lines 51-52: Pulse infrastructure; implemented by `H=hamiltonian(assume(spin_system,'nmr'))`.
- Lines 58-59: WURST chirp waveform in amplitude-frequency coordinates; implemented by `[~,~,durs,~,amps,~,frqs]=chirp_pulse(100,0.1,2000,16,'wurst')`.
- Lines 61-62: Chirp frequency shift; implemented by `frqs=frqs+1000`.
- Lines 64-66: Soft pulse in amplitude-frequency coordinates; implemented by `parameters.rho0=shaped_pulse_af(spin_system,H,Lx,Ly,parameters.rho0, frqs,amps,durs,pi/2,2)`.
- Lines 68-69: Homospoil gradient to kill stray transverse magnetisation; implemented by `parameters.rho0=homospoil(spin_system,parameters.rho0,'destroy')`.
- Lines 71-72: Global hard pulse to detect longitudinal magnetisation; implemented by `parameters.rho0=step(spin_system,Ly,parameters.rho0,pi/2)`.
- Lines 74-75: Acquisition; implemented by `fid=acquire(spin_system,parameters,H,R,K)`.
- Lines 77-78: Apodisation; implemented by `fid=apodisation(spin_system,fid,{{'exp',6}})`.
- Lines 80-81: Fourier transform; implemented by `spectrum=fftshift(fft(fid,parameters.zerofill))`.
- Lines 83-84: Plotting; implemented by `kfigure(); plot_1d(spin_system,real(spectrum),parameters)`.

### Control flow inferred from the code

- Line 26: `for` loop over `n=1:30`.

### Key state/data transformations

- Lines 14: computes `sys.magnet` using `sys.magnet=14.1`.
- Lines 17-19: computes `sys.isotopes` using `sys.isotopes={'1H','1H','1H','1H','1H','1H','1H','1H','1H','1H','1H', '1H','1H','1H','1H','1H','1H','1H','1H','1H','1H','1H', '1H','1H','1H','1H','1H','1H','1H','1H','1H…`.
- Lines 22: computes `inter.zeeman.scalar` using `inter.zeeman.scalar=num2cell(linspace(-4,4,31))`.
- Lines 25: computes `inter.coupling.scalar` using `inter.coupling.scalar=cell(31)`.
- Lines 27: computes `inter.coupling.scalar{n,n+1}` using `inter.coupling.scalar{n,n+1}=20`.
- Lines 31: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 32: computes `bas.approximation` using `bas.approximation='IK-2'`.
- Lines 33: computes `bas.connectivity` using `bas.connectivity='scalar_couplings'`.
- Lines 34: computes `bas.space_level` using `bas.space_level=1`.
- Lines 37: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 41: computes `parameters.spins` using `parameters.spins={'1H'}`.
- Lines 42: computes `parameters.rho0` using `parameters.rho0=state(spin_system,'Lz','1H')`.
- Lines 43: computes `parameters.coil` using `parameters.coil=state(spin_system,'L+','1H')`.
- Lines 44: computes `parameters.decouple` using `parameters.decouple={}`.
- Lines 45: computes `parameters.offset` using `parameters.offset=0`.
- Lines 46: computes `parameters.sweep` using `parameters.sweep=5100`.
- Lines 47: computes `parameters.npoints` using `parameters.npoints=2048`.
- Lines 48: computes `parameters.zerofill` using `parameters.zerofill=8192`.

## Implementation structure

- Chirp inversion pulse using the Fokker-Planck formalism. Fewer
- points are required by the amplitude-frequency method than the
- "two points per period of the largest frequency" Nyquist-Shan-
- non condition would need for the {Cx,Cy} parameterised simula-
- tion of a chirped pulse.
- Calculation time: seconds
- Magnetic field
- Isotopes
- Zeeman interactions
- Couplings
- Basis set
- Spinach housekeeping

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `num2cell()`, `create()`, `basis()`, `state()`, `hamiltonian()`, `assume()`, `relaxation()`, `kinetics()`, `operator()`, `chirp_pulse()`, `shaped_pulse_af()`, `homospoil()`, `step()`, `acquire()`, `apodisation()`, `fftshift()`.
