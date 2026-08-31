# examples/shaped_pulses/shaped_pulse_fp.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/shaped_pulses/shaped_pulse_fp.m`
- Signature: `shaped_pulse_fp()`
- Total lines: 79

## Purpose

An off-resonance rectangular soft pulse simulated using the Fokker-Planck formalism. Note that the pulse frequency off- set accumulates as additional phase during the pulse in the same way as it would during a chirp. Calculation time: seconds

## Physical / mathematical content

- Shaped-pulse examples. These scripts demonstrate amplitude, phase, frequency, and gradient waveform design, including adiabatic sweeps, excitation profiles, and hardware-response considerations.
- The file uses a Fokker-Planck-style enlarged state space in which spatial or orientational coordinates are promoted to extra dimensions and coupled to spin dynamics through differential operators.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 12-13: Magnetic field; implemented by `sys.magnet=14.1`.
- Lines 15-18: Isotopes; implemented by `sys.isotopes={'1H','1H','1H','1H','1H','1H','1H','1H','1H','1H','1H', '1H','1H','1H','1H','1H','1H','1H','1H','1H','1H','1H', '1H','1H','1H','1H','1H','1H','1H','1H','1H…`.
- Lines 20-21: Zeeman interactions; implemented by `inter.zeeman.scalar=num2cell(linspace(-4,4,31))`.
- Lines 23-24: Couplings; implemented by `inter.coupling.scalar=cell(31)`.
- Lines 29-30: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 35-36: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 39-40: Background Hamiltonian; implemented by `H=hamiltonian(assume(spin_system,'nmr'))`.
- Lines 42-43: Control operators; implemented by `Lx=operator(spin_system,'Lx','1H')`.
- Lines 46-47: Initial condition; implemented by `rho=state(spin_system,'Lz','1H')`.
- Lines 49-50: Fokker-Planck soft pulse at 1922.4 Hz frequency offset; implemented by `rho=shaped_pulse_af(spin_system,H,Lx,Ly,rho,1922.4,50.0,5e-3,-pi/2,2)`.
- Lines 52-53: Set up acquisition; implemented by `parameters.spins={'1H'}`.
- Lines 63-64: Run acquisition; implemented by `fid=liquid(spin_system,@acquire,parameters,'nmr')`.
- Lines 66-67: Phasing; implemented by `fid=fid*exp(-1i*0.67)`.
- Lines 69-70: Apodisation; implemented by `fid=apodisation(spin_system,fid,{{'exp',6}})`.
- Lines 72-73: Fourier transform; implemented by `spectrum=fftshift(fft(fid,parameters.zerofill))`.
- Lines 75-76: Plotting; implemented by `kfigure(); plot_1d(spin_system,real(spectrum),parameters)`.

### Control flow inferred from the code

- Line 25: `for` loop over `n=1:30`.

### Key state/data transformations

- Lines 13: computes `sys.magnet` using `sys.magnet=14.1`.
- Lines 16-18: computes `sys.isotopes` using `sys.isotopes={'1H','1H','1H','1H','1H','1H','1H','1H','1H','1H','1H', '1H','1H','1H','1H','1H','1H','1H','1H','1H','1H','1H', '1H','1H','1H','1H','1H','1H','1H','1H','1H…`.
- Lines 21: computes `inter.zeeman.scalar` using `inter.zeeman.scalar=num2cell(linspace(-4,4,31))`.
- Lines 24: computes `inter.coupling.scalar` using `inter.coupling.scalar=cell(31)`.
- Lines 26: computes `inter.coupling.scalar{n,n+1}` using `inter.coupling.scalar{n,n+1}=10`.
- Lines 30: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 31: computes `bas.approximation` using `bas.approximation='IK-2'`.
- Lines 32: computes `bas.connectivity` using `bas.connectivity='scalar_couplings'`.
- Lines 33: computes `bas.space_level` using `bas.space_level=1`.
- Lines 36: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 40: computes `H` using `H=hamiltonian(assume(spin_system,'nmr'))`.
- Lines 43: computes `Lx` using `Lx=operator(spin_system,'Lx','1H')`.
- Lines 44: computes `Ly` using `Ly=operator(spin_system,'Ly','1H')`.
- Lines 47: computes `rho` using `rho=state(spin_system,'Lz','1H')`.
- Lines 53: computes `parameters.spins` using `parameters.spins={'1H'}`.
- Lines 54: computes `parameters.rho0` using `parameters.rho0=rho`.
- Lines 55: computes `parameters.coil` using `parameters.coil=state(spin_system,'L+','1H')`.
- Lines 56: computes `parameters.decouple` using `parameters.decouple={}`.

## Implementation structure

- An off-resonance rectangular soft pulse simulated using the
- Fokker-Planck formalism. Note that the pulse frequency off-
- set accumulates as additional phase during the pulse in the
- same way as it would during a chirp.
- Calculation time: seconds
- Magnetic field
- Isotopes
- Zeeman interactions
- Couplings
- Basis set
- Spinach housekeeping
- Background Hamiltonian

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `num2cell()`, `create()`, `basis()`, `hamiltonian()`, `assume()`, `operator()`, `state()`, `shaped_pulse_af()`, `liquid()`, `apodisation()`, `fftshift()`, `kfigure()`, `plot_1d()`.
