# examples/shaped_pulses/shaped_pulse_vg.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/shaped_pulses/shaped_pulse_vg.m`
- Signature: `shaped_pulse_vg()`
- Total lines: 83

## Purpose

Veshtort-Griffin E1000B 90-degree selective pulse applied to a system of 31 proton spins with nearest-neighbor J co- uplings and linear coupling topology. Calculation time: seconds

## Physical / mathematical content

- Shaped-pulse examples. These scripts demonstrate amplitude, phase, frequency, and gradient waveform design, including adiabatic sweeps, excitation profiles, and hardware-response considerations.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 11-12: Magnetic field; implemented by `sys.magnet=14.1`.
- Lines 14-17: Isotopes; implemented by `sys.isotopes={ '1H','1H','1H','1H','1H','1H','1H','1H','1H','1H','1H', '1H','1H','1H','1H','1H','1H','1H','1H','1H','1H','1H', '1H','1H','1H','1H','1H','1H','1H','1H','1…`.
- Lines 19-20: Zeeman interactions; implemented by `inter.zeeman.scalar=num2cell(linspace(-4,4,31))`.
- Lines 22-23: Couplings; implemented by `inter.coupling.scalar=cell(31)`.
- Lines 28-29: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 34-35: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 38-39: Assumptions; implemented by `spin_system=assume(spin_system,'nmr')`.
- Lines 41-42: Hamiltonian superoperator; implemented by `H=hamiltonian(spin_system)`.
- Lines 44-45: Control and offset operators; implemented by `Lx=operator(spin_system,'Lx','1H')`.
- Lines 48-49: Initial state; implemented by `rho=state(spin_system,'Lz','1H')`.
- Lines 51-52: Pulse parameters; implemented by `duration=1e-2; time_grid=duration*ones(1,500)/500`.
- Lines 55-57: Pulse execution with 480 Hz frequency offset; implemented by `rho=shaped_pulse_xy(spin_system,H+2*pi*480*Lz, {Lx},{amplitudes},time_grid,rho,'expv-pwc')`.
- Lines 59-60: Set up acquisition; implemented by `parameters.spins={'1H'}`.
- Lines 70-71: Run acquisition; implemented by `fid=liquid(spin_system,@acquire,parameters,'nmr')`.
- Lines 73-74: Apodisation; implemented by `fid=apodisation(spin_system,fid,{{'exp',6}})`.
- Lines 76-77: Fourier transform; implemented by `spectrum=fftshift(fft(fid,parameters.zerofill))`.
- Lines 79-80: Plotting; implemented by `kfigure(); plot_1d(spin_system,imag(spectrum),parameters)`.

### Control flow inferred from the code

- Line 24: `for` loop over `n=1:30`.

### Key state/data transformations

- Lines 12: computes `sys.magnet` using `sys.magnet=14.1`.
- Lines 15-17: computes `sys.isotopes` using `sys.isotopes={ '1H','1H','1H','1H','1H','1H','1H','1H','1H','1H','1H', '1H','1H','1H','1H','1H','1H','1H','1H','1H','1H','1H', '1H','1H','1H','1H','1H','1H','1H','1H','1…`.
- Lines 20: computes `inter.zeeman.scalar` using `inter.zeeman.scalar=num2cell(linspace(-4,4,31))`.
- Lines 23: computes `inter.coupling.scalar` using `inter.coupling.scalar=cell(31)`.
- Lines 25: computes `inter.coupling.scalar{n,n+1}` using `inter.coupling.scalar{n,n+1}=10`.
- Lines 29: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 30: computes `bas.approximation` using `bas.approximation='IK-2'`.
- Lines 31: computes `bas.connectivity` using `bas.connectivity='scalar_couplings'`.
- Lines 32: computes `bas.space_level` using `bas.space_level=1`.
- Lines 35: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 42: computes `H` using `H=hamiltonian(spin_system)`.
- Lines 45: computes `Lx` using `Lx=operator(spin_system,'Lx','1H')`.
- Lines 46: computes `Lz` using `Lz=operator(spin_system,'Lz','1H')`.
- Lines 49: computes `rho` using `rho=state(spin_system,'Lz','1H')`.
- Lines 52: computes `duration` using `duration=1e-2; time_grid=duration*ones(1,500)/500`.
- Lines 53: computes `amplitudes` using `amplitudes=vg_pulse('E1000B',500,duration)`.
- Lines 60: computes `parameters.spins` using `parameters.spins={'1H'}`.
- Lines 61: computes `parameters.rho0` using `parameters.rho0=rho`.

## Implementation structure

- Veshtort-Griffin E1000B 90-degree selective pulse applied
- to a system of 31 proton spins with nearest-neighbor J co-
- uplings and linear coupling topology.
- Calculation time: seconds
- Magnetic field
- Isotopes
- Zeeman interactions
- Couplings
- Basis set
- Spinach housekeeping
- Assumptions
- Hamiltonian superoperator

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `num2cell()`, `create()`, `basis()`, `assume()`, `hamiltonian()`, `operator()`, `state()`, `vg_pulse()`, `shaped_pulse_xy()`, `liquid()`, `apodisation()`, `fftshift()`, `kfigure()`, `plot_1d()`.
