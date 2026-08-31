# examples/nmr_spen/dosy_oneshot_1.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_spen/dosy_oneshot_1.m`
- Signature: `dosy_oneshot_1()`
- Total lines: 88

## Purpose

Oneshot DOSY pulse sequence for a system of three coupled spins with different relaxation rates. Timing: minutes on NVidia Tesla A100, much longer on CPU

## Physical / mathematical content

- SPEN / ultrafast NMR examples. These files encode spatially dependent phase evolution and acquisition, linking pulse gradients, diffusion attenuation, and single-scan multidimensional encoding.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 10-11: Magnetic field; implemented by `sys.magnet=11.7428`.
- Lines 13-14: Spin system; implemented by `sys.isotopes={'1H','1H','1H'}`.
- Lines 21-22: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 25-26: Relaxation theory parameters; implemented by `inter.relaxation={'t1_t2'}`.
- Lines 32-33: Algorithmic options; implemented by `sys.disable={'pt'}`.
- Lines 36-37: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 40-41: Acquisition parameters; implemented by `parameters.spins={'1H'}`.
- Lines 48-49: Sample geometry; implemented by `parameters.dims=0.015`.
- Lines 53-54: Relaxation phantom; implemented by `parameters.rlx_ph={ones(parameters.npts,1)}`.
- Lines 57-58: Initial and detection state phantoms; implemented by `parameters.rho0_ph={ones(parameters.npts,1)}`.
- Lines 63-64: Diffusion; implemented by `parameters.diff=18.55e-10`.
- Lines 66-67: Gradient parameters; implemented by `parameters.g_amp=0.255`.
- Lines 70-71: Timing parameters; implemented by `parameters.g_dur=0.001`.
- Lines 75-76: Call the sequence in the imaging context; implemented by `fid=imaging(spin_system,@dosy_oneshot,parameters)`.
- Lines 78-79: Apodisation; implemented by `fid=apodisation(spin_system,fid,{{'exp',5}})`.
- Lines 81-82: Fourier transform; implemented by `spectrum=fftshift(fft(fid,parameters.zerofill))`.
- Lines 84-85: Plotting; implemented by `kfigure(); plot_1d(spin_system,-real(spectrum),parameters)`.

### Key state/data transformations

- Lines 11: computes `sys.magnet` using `sys.magnet=11.7428`.
- Lines 14: computes `sys.isotopes` using `sys.isotopes={'1H','1H','1H'}`.
- Lines 15: computes `inter.zeeman.scalar` using `inter.zeeman.scalar={4.70 3.50 1.50}`.
- Lines 16: computes `inter.coupling.scalar{1,2}` using `inter.coupling.scalar{1,2}=15`.
- Lines 17: computes `inter.coupling.scalar{2,3}` using `inter.coupling.scalar{2,3}=15`.
- Lines 18: computes `inter.coupling.scalar{1,3}` using `inter.coupling.scalar{1,3}=10`.
- Lines 19: computes `inter.coupling.scalar{3,3}` using `inter.coupling.scalar{3,3}=0`.
- Lines 22: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 23: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 26: computes `inter.relaxation` using `inter.relaxation={'t1_t2'}`.
- Lines 27: computes `inter.rlx_keep` using `inter.rlx_keep='diagonal'`.
- Lines 28: computes `inter.equilibrium` using `inter.equilibrium='zero'`.
- Lines 29: computes `inter.r1_rates` using `inter.r1_rates=num2cell(1./[0.1952 0.2100 0.2500])`.
- Lines 30: computes `inter.r2_rates` using `inter.r2_rates=num2cell(1./[0.1602 0.1802 0.1902])`.
- Lines 33: computes `sys.disable` using `sys.disable={'pt'}`.
- Lines 34: computes `sys.enable` using `sys.enable={'greedy'}`.
- Lines 37: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 41: computes `parameters.spins` using `parameters.spins={'1H'}`.

## Implementation structure

- Oneshot DOSY pulse sequence for a system of three coupled
- spins with different relaxation rates.
- Timing: minutes on NVidia Tesla A100, much longer on CPU
- Magnetic field
- Spin system
- Basis set
- Relaxation theory parameters
- Algorithmic options
- Spinach housekeeping
- Acquisition parameters
- Sample geometry
- Relaxation phantom

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `num2cell()`, `create()`, `basis()`, `relaxation()`, `state()`, `imaging()`, `apodisation()`, `fftshift()`, `kfigure()`, `plot_1d()`.
