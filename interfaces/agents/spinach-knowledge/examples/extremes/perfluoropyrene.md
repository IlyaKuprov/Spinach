# examples/extremes/perfluoropyrene.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/extremes/perfluoropyrene.m`
- Signature: `perfluoropyrene()`
- Total lines: 64

## Purpose

X-band pulsed ESR spectrum of perfluoropyrene cation radical, computed using brute force operator algebra in the full 4,194,304 -dimensional Liouville space. This is deliberate -a much faster calculation is, of course, possible with a restricted basis set. This calculation requires at least 64GB of RAM and illustrates the per- formance of trajectory-level state space restriction in Spinach. Calculation time: minutes

## Physical / mathematical content

- Extreme-regime examples. These scripts exercise Spinach in unusually large, stiff, high-field, low-field, or otherwise numerically demanding regimes where approximations, conditioning, and basis-size control are central.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 15-16: Ignore coordinate information (HFCs provided); implemented by `options.no_xyz=1`.
- Lines 18-20: Read the spin system properties (vacuum DFT calculation); implemented by `[sys,inter]=g2spinach(gparse('../standard_systems/perfluoropyrene_cation.log'), {{'E','E'},{'F','19F'}},[0 0],options)`.
- Lines 21-22: Magnet induction; implemented by `sys.magnet=0.33`.
- Lines 24-25: Relaxation theory; implemented by `inter.relaxation={'damp'}`.
- Lines 30-31: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 34-35: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 38-39: Set the sequence parameters; implemented by `parameters.spins={'E'}`.
- Lines 51-52: Simulation; implemented by `fid=liquid(spin_system,@acquire,parameters,'esr')`.
- Lines 54-55: Apodisation; implemented by `fid=apodisation(spin_system,fid,{{'none'}})`.
- Lines 57-58: Fourier transform; implemented by `spectrum=fftshift(fft(fid,parameters.zerofill))`.
- Lines 60-61: Plotting; implemented by `kfigure(); plot_1d(spin_system,real(spectrum),parameters)`.

### Key state/data transformations

- Lines 16: computes `options.no_xyz` using `options.no_xyz=1`.
- Lines 19-20: computes `[sys,inter]` using `[sys,inter]=g2spinach(gparse('../standard_systems/perfluoropyrene_cation.log'), {{'E','E'},{'F','19F'}},[0 0],options)`.
- Lines 22: computes `sys.magnet` using `sys.magnet=0.33`.
- Lines 25: computes `inter.relaxation` using `inter.relaxation={'damp'}`.
- Lines 26: computes `inter.rlx_keep` using `inter.rlx_keep='diagonal'`.
- Lines 27: computes `inter.equilibrium` using `inter.equilibrium='zero'`.
- Lines 28: computes `inter.damp_rate` using `inter.damp_rate=2e6`.
- Lines 31: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 32: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 35: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 39: computes `parameters.spins` using `parameters.spins={'E'}`.
- Lines 40: computes `parameters.rho0` using `parameters.rho0=state(spin_system,'L+','E')`.
- Lines 41: computes `parameters.coil` using `parameters.coil=state(spin_system,'L+','E')`.
- Lines 42: computes `parameters.decouple` using `parameters.decouple={}`.
- Lines 43: computes `parameters.offset` using `parameters.offset=0`.
- Lines 44: computes `parameters.sweep` using `parameters.sweep=3e8`.
- Lines 45: computes `parameters.npoints` using `parameters.npoints=2048`.
- Lines 46: computes `parameters.zerofill` using `parameters.zerofill=4096`.

## Implementation structure

- X-band pulsed ESR spectrum of perfluoropyrene cation radical, computed
- using brute force operator algebra in the full 4,194,304 -dimensional
- Liouville space. This is deliberate -a much faster calculation is, of
- course, possible with a restricted basis set.
- This calculation requires at least 64GB of RAM and illustrates the per-
- formance of trajectory-level state space restriction in Spinach.
- Calculation time: minutes
- Ignore coordinate information (HFCs provided)
- Read the spin system properties (vacuum DFT calculation)
- Magnet induction
- Relaxation theory
- Basis set

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `g2spinach()`, `gparse()`, `create()`, `basis()`, `state()`, `liquid()`, `apodisation()`, `fftshift()`, `kfigure()`, `plot_1d()`.
