# examples/shaped_pulses/shaped_pulse_slr.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/shaped_pulses/shaped_pulse_slr.m`
- Signature: `shaped_pulse_slr()`
- Total lines: 81

## Purpose

Shinnar-Le Roux band-selective 90-degree pulse on a chain of 31 strongly coupled protons. Calculation time: seconds

## Physical / mathematical content

- Shaped-pulse examples. These scripts demonstrate amplitude, phase, frequency, and gradient waveform design, including adiabatic sweeps, excitation profiles, and hardware-response considerations.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 10-11: Magnetic field; implemented by `sys.magnet=14.1`.
- Lines 13-14: Isotopes; implemented by `sys.isotopes=repmat({'1H'},1,31)`.
- Lines 16-17: Zeeman interactions; implemented by `inter.zeeman.scalar=num2cell(linspace(-4,4,31))`.
- Lines 19-20: Couplings; implemented by `inter.coupling.scalar=cell(31)`.
- Lines 25-26: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 31-32: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 36-37: Hamiltonian and control operators; implemented by `H=hamiltonian(spin_system)`.
- Lines 41-42: Initial state; implemented by `rho=state(spin_system,'Lz','1H')`.
- Lines 44-45: Inversion SLR pulse waveform; implemented by `[Cx,Cy,durs]=slr_pulse(256,15e-3,32,pi/2,0.01,0.01)`.
- Lines 47-49: Run the SLR pulse; implemented by `rho=shaped_pulse_xy(spin_system,H,{Lx,Ly},{Cx,Cy}, durs,rho,'expv-pwc')`.
- Lines 51-52: Set up acquisition; implemented by `parameters.spins={'1H'}`.
- Lines 62-63: Run acquisition; implemented by `fid=liquid(spin_system,@acquire,parameters,'nmr')`.
- Lines 65-66: Apodisation; implemented by `fid=apodisation(spin_system,fid,{{'exp',6}})`.
- Lines 68-69: Fourier transform; implemented by `spectrum=fftshift(fft(fid,parameters.zerofill))`.
- Lines 71-72: Plotting; implemented by `kfigure(); scale_figure([2.0 1.0])`.

### Control flow inferred from the code

- Line 21: `for` loop over `n=1:30`.

### Key state/data transformations

- Lines 11: computes `sys.magnet` using `sys.magnet=14.1`.
- Lines 14: computes `sys.isotopes` using `sys.isotopes=repmat({'1H'},1,31)`.
- Lines 17: computes `inter.zeeman.scalar` using `inter.zeeman.scalar=num2cell(linspace(-4,4,31))`.
- Lines 20: computes `inter.coupling.scalar` using `inter.coupling.scalar=cell(31)`.
- Lines 22: computes `inter.coupling.scalar{n,n+1}` using `inter.coupling.scalar{n,n+1}=10`.
- Lines 26: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 27: computes `bas.approximation` using `bas.approximation='IK-2'`.
- Lines 28: computes `bas.connectivity` using `bas.connectivity='scalar_couplings'`.
- Lines 29: computes `bas.space_level` using `bas.space_level=1`.
- Lines 32: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 37: computes `H` using `H=hamiltonian(spin_system)`.
- Lines 38: computes `Lx` using `Lx=operator(spin_system,'Lx','1H')`.
- Lines 39: computes `Ly` using `Ly=operator(spin_system,'Ly','1H')`.
- Lines 42: computes `rho` using `rho=state(spin_system,'Lz','1H')`.
- Lines 45: computes `[Cx,Cy,durs]` using `[Cx,Cy,durs]=slr_pulse(256,15e-3,32,pi/2,0.01,0.01)`.
- Lines 52: computes `parameters.spins` using `parameters.spins={'1H'}`.
- Lines 53: computes `parameters.rho0` using `parameters.rho0=rho`.
- Lines 54: computes `parameters.coil` using `parameters.coil=state(spin_system,'L+','1H')`.

## Implementation structure

- Shinnar-Le Roux band-selective 90-degree pulse on a chain of
- 31 strongly coupled protons.
- Calculation time: seconds
- Magnetic field
- Isotopes
- Zeeman interactions
- Couplings
- Basis set
- Spinach housekeeping
- Hamiltonian and control operators
- Initial state
- Inversion SLR pulse waveform

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `num2cell()`, `create()`, `basis()`, `assume()`, `hamiltonian()`, `operator()`, `state()`, `slr_pulse()`, `shaped_pulse_xy()`, `liquid()`, `apodisation()`, `fftshift()`, `kfigure()`, `scale_figure()`, `subplot()`, `cumsum()`.
