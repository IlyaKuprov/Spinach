# examples/shaped_pulses/shaped_pulse_q5.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/shaped_pulses/shaped_pulse_q5.m`
- Signature: `shaped_pulse_q5()`
- Total lines: 91

## Purpose

90-degree Q5 pulse on a chain of 31 strongly coupled protons. Calculation time: seconds

## Physical / mathematical content

- Shaped-pulse examples. These scripts demonstrate amplitude, phase, frequency, and gradient waveform design, including adiabatic sweeps, excitation profiles, and hardware-response considerations.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 9-10: Magnetic field; implemented by `sys.magnet=14.1`.
- Lines 12-15: Isotopes; implemented by `sys.isotopes={'1H','1H','1H','1H','1H','1H','1H','1H','1H','1H','1H', '1H','1H','1H','1H','1H','1H','1H','1H','1H','1H','1H', '1H','1H','1H','1H','1H','1H','1H','1H','1H…`.
- Lines 17-18: Zeeman interactions; implemented by `inter.zeeman.scalar=num2cell(linspace(-4,4,31))`.
- Lines 20-21: Couplings; implemented by `inter.coupling.scalar=cell(31)`.
- Lines 26-27: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 32-33: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 36-37: Assumptions; implemented by `spin_system=assume(spin_system,'nmr')`.
- Lines 39-40: Hamiltonian superoperator; implemented by `H=hamiltonian(spin_system)`.
- Lines 42-43: Control and offset operators; implemented by `Lx=operator(spin_system,'Lx','1H')`.
- Lines 47-48: Initial state; implemented by `rho=state(spin_system,'Lz','1H')`.
- Lines 50-51: Q5 pulse setup; implemented by `npoints=200; duration=0.012`.
- Lines 55-56: Pulse calibration; implemented by `integral=sum(A)`.
- Lines 60-61: Convert into Cartesian coordinates; implemented by `[Cx,Cy]=polar2cartesian(A,phi)`.
- Lines 63-65: Pulse execution with a 480 Hz frequency offset; implemented by `rho=shaped_pulse_xy(spin_system,H+2*pi*480*Lz, {Lx,Ly},{Cx,Cy},time_grid,rho,'expv-pwc')`.
- Lines 67-68: Set up acquisition; implemented by `parameters.spins={'1H'}`.
- Lines 78-79: Run acquisition; implemented by `fid=liquid(spin_system,@acquire,parameters,'nmr')`.
- Lines 81-82: Apodisation; implemented by `fid=apodisation(spin_system,fid,{{'exp',6}})`.
- Lines 84-85: Fourier transform; implemented by `spectrum=fftshift(fft(fid,parameters.zerofill))`.

### Control flow inferred from the code

- Line 22: `for` loop over `n=1:30`.

### Key state/data transformations

- Lines 10: computes `sys.magnet` using `sys.magnet=14.1`.
- Lines 13-15: computes `sys.isotopes` using `sys.isotopes={'1H','1H','1H','1H','1H','1H','1H','1H','1H','1H','1H', '1H','1H','1H','1H','1H','1H','1H','1H','1H','1H','1H', '1H','1H','1H','1H','1H','1H','1H','1H','1H…`.
- Lines 18: computes `inter.zeeman.scalar` using `inter.zeeman.scalar=num2cell(linspace(-4,4,31))`.
- Lines 21: computes `inter.coupling.scalar` using `inter.coupling.scalar=cell(31)`.
- Lines 23: computes `inter.coupling.scalar{n,n+1}` using `inter.coupling.scalar{n,n+1}=10`.
- Lines 27: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 28: computes `bas.approximation` using `bas.approximation='IK-2'`.
- Lines 29: computes `bas.connectivity` using `bas.connectivity='scalar_couplings'`.
- Lines 30: computes `bas.space_level` using `bas.space_level=1`.
- Lines 33: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 40: computes `H` using `H=hamiltonian(spin_system)`.
- Lines 43: computes `Lx` using `Lx=operator(spin_system,'Lx','1H')`.
- Lines 44: computes `Ly` using `Ly=operator(spin_system,'Ly','1H')`.
- Lines 45: computes `Lz` using `Lz=operator(spin_system,'Lz','1H')`.
- Lines 48: computes `rho` using `rho=state(spin_system,'Lz','1H')`.
- Lines 51: computes `npoints` using `npoints=200; duration=0.012`.
- Lines 52: computes `time_grid` using `time_grid=duration*ones(1,npoints)/npoints`.
- Lines 53: computes `[A,phi]` using `[A,phi]=read_wave('q5_1000.pk',npoints)`.

## Implementation structure

- 90-degree Q5 pulse on a chain of 31 strongly coupled protons.
- Calculation time: seconds
- Magnetic field
- Isotopes
- Zeeman interactions
- Couplings
- Basis set
- Spinach housekeeping
- Assumptions
- Hamiltonian superoperator
- Control and offset operators
- Initial state

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `num2cell()`, `create()`, `basis()`, `assume()`, `hamiltonian()`, `operator()`, `state()`, `read_wave()`, `polar2cartesian()`, `shaped_pulse_xy()`, `liquid()`, `apodisation()`, `fftshift()`, `kfigure()`, `plot_1d()`.
