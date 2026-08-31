# examples/nmr_spen/ufmq_2spin.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_spen/ufmq_2spin.m`
- Signature: `ufmq_2spin()`
- Total lines: 106

## Purpose

2Q ultrafast MaxQ NMR spectrum for a coupled two-spin system in the presence of realistic diffusion. Calculation time: minutes on NVidia Tesla A100, much longer on CPU

## Physical / mathematical content

- SPEN / ultrafast NMR examples. These files encode spatially dependent phase evolution and acquisition, linking pulse gradients, diffusion attenuation, and single-scan multidimensional encoding.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 13-14: Magnetic field; implemented by `sys.magnet=14.1`.
- Lines 16-17: Chemical shifts; implemented by `sys.isotopes={'1H','1H'}`.
- Lines 20-21: J-coupling; implemented by `inter.coupling.scalar{1,2}=8.0`.
- Lines 24-25: Coherence selection; implemented by `parameters.mqorder=+2`.
- Lines 27-28: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 31-32: Algorithmic options; implemented by `sys.disable={'pt'}`.
- Lines 35-36: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 39-40: Sample geometry; implemented by `parameters.dims=0.015`.
- Lines 44-45: Relaxation phantom; implemented by `parameters.rlx_ph={}`.
- Lines 48-49: Initial and detection state phantoms; implemented by `parameters.rho0_ph={ones(parameters.npts,1)}`.
- Lines 54-55: Diffusion and flow; implemented by `parameters.u=zeros(parameters.npts,1)`.
- Lines 58-59: Acquisition parameters; implemented by `parameters.spins={'1H'}`.
- Lines 65-66: Calculate maximal k-value; implemented by `k_max=parameters.npoints/parameters.dims`.
- Lines 68-69: Get acquisition gradient duration; implemented by `Ta=parameters.deltat*parameters.npoints`.
- Lines 71-72: Get acquisition gradient amplitude; implemented by `parameters.Ga=((2*pi)*k_max)/(spin(parameters.spins{1})*Ta)`.
- Lines 74-75: Timings; implemented by `parameters.delay=0.041`.
- Lines 77-78: Encoding parameters; implemented by `parameters.pulsenpoints=500`.
- Lines 85-86: Simulation; implemented by `ktdata=imaging(spin_system,@ufmq,parameters)`.

### Key state/data transformations

- Lines 14: computes `sys.magnet` using `sys.magnet=14.1`.
- Lines 17: computes `sys.isotopes` using `sys.isotopes={'1H','1H'}`.
- Lines 18: computes `inter.zeeman.scalar` using `inter.zeeman.scalar={0.50,0.15}`.
- Lines 21: computes `inter.coupling.scalar{1,2}` using `inter.coupling.scalar{1,2}=8.0`.
- Lines 22: computes `inter.coupling.scalar{2,2}` using `inter.coupling.scalar{2,2}=0`.
- Lines 25: computes `parameters.mqorder` using `parameters.mqorder=+2`.
- Lines 28: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 29: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 32: computes `sys.disable` using `sys.disable={'pt'}`.
- Lines 33: computes `sys.enable` using `sys.enable={'greedy'}`.
- Lines 36: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 40: computes `parameters.dims` using `parameters.dims=0.015`.
- Lines 41: computes `parameters.npts` using `parameters.npts=500`.
- Lines 42: computes `parameters.deriv` using `parameters.deriv={'period',7}`.
- Lines 45: computes `parameters.rlx_ph` using `parameters.rlx_ph={}`.
- Lines 46: computes `parameters.rlx_op` using `parameters.rlx_op={}`.
- Lines 49: computes `parameters.rho0_ph` using `parameters.rho0_ph={ones(parameters.npts,1)}`.
- Lines 50: computes `parameters.rho0_st` using `parameters.rho0_st={state(spin_system,'Lz','1H')}`.

## Implementation structure

- 2Q ultrafast MaxQ NMR spectrum for a coupled two-spin
- system in the presence of realistic diffusion.
- Calculation time: minutes on NVidia Tesla A100,
- much longer on CPU
- Magnetic field
- Chemical shifts
- J-coupling
- Coherence selection
- Basis set
- Algorithmic options
- Spinach housekeeping
- Sample geometry

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `state()`, `spin()`, `imaging()`, `kfigure()`, `subplot()`, `fftshift()`, `plot_uf()`.
