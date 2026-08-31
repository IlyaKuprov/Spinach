# examples/nmr_liquids/crazed_test.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_liquids/crazed_test.m`
- Signature: `crazed_test()`
- Total lines: 60

## Purpose

Long range intermolecular coherences predicted by Warren and co-workers Calculation time: seconds

## Physical / mathematical content

- Liquid-state NMR examples. The physics is scalar-coupling-mediated coherence transfer in weakly or moderately coupled spin systems, often in Liouville space. Typical mechanisms include INEPT-style polarisation transfer, J-refocusing, phase cycling, indirect evolution, and multidimensional detection.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 13-14: Specify system parameters; implemented by `sys.magnet=6.0`.
- Lines 23-24: Use the complete basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 27-28: Spinach code; implemented by `spin_system=create(sys,inter)`.
- Lines 31-32: Sequence parameters; implemented by `parameters.angle=pi/2`.
- Lines 40-41: Thermal equilibrium state; implemented by `[H,Q]=hamiltonian(assume(spin_system,'labframe'),'left')`.
- Lines 44-45: CRAZED simulation; implemented by `fid=crystal(spin_system,@crazed,parameters,'nmr')`.
- Lines 47-48: Apodisation; implemented by `fid=apodisation(spin_system,fid,{{'cos'},{'cos'}})`.
- Lines 50-52: Fourier transform; implemented by `spectrum=fftshift(fft2(fid,parameters.zerofill(2), parameters.zerofill(1)))`.
- Lines 54-55: Plotting; implemented by `kfigure(); scale_figure([1.5 2.0])`.

### Key state/data transformations

- Lines 14: computes `sys.magnet` using `sys.magnet=6.0`.
- Lines 15: computes `sys.isotopes` using `sys.isotopes={'1H','1H','1H','1H'}`.
- Lines 16: computes `inter.zeeman.scalar` using `inter.zeeman.scalar={2.0,2.0,8.0,8.0}`.
- Lines 17: computes `inter.coordinates` using `inter.coordinates={[0 0 0]`.
- Lines 21: computes `inter.temperature` using `inter.temperature=100`.
- Lines 24: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 25: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 28: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 32: computes `parameters.angle` using `parameters.angle=pi/2`.
- Lines 33: computes `parameters.offset` using `parameters.offset=1300`.
- Lines 34: computes `parameters.sweep` using `parameters.sweep=5000`.
- Lines 35: computes `parameters.npoints` using `parameters.npoints=[512 512]`.
- Lines 36: computes `parameters.zerofill` using `parameters.zerofill=[2048 2048]`.
- Lines 37: computes `parameters.spins` using `parameters.spins={'1H'}`.
- Lines 38: computes `parameters.orientation` using `parameters.orientation=[0 0 0]`.
- Lines 41: computes `[H,Q]` using `[H,Q]=hamiltonian(assume(spin_system,'labframe'),'left')`.
- Lines 42: computes `parameters.rho0` using `parameters.rho0=equilibrium(spin_system,H,Q,[0 0 0])`.
- Lines 45: computes `fid` using `fid=crystal(spin_system,@crazed,parameters,'nmr')`.

## Implementation structure

- Long range intermolecular coherences predicted by Warren and
- co-workers
- Calculation time: seconds
- Specify system parameters
- Use the complete basis set
- Spinach code
- Sequence parameters
- Thermal equilibrium state
- CRAZED simulation
- Apodisation
- Fourier transform
- Plotting

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `hamiltonian()`, `assume()`, `equilibrium()`, `crystal()`, `apodisation()`, `fftshift()`, `fft2()`, `kfigure()`, `scale_figure()`, `plot_2d()`.
