# examples/nmr_spen/ufmq_4spin.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_spen/ufmq_4spin.m`
- Signature: `ufmq_4spin()`
- Total lines: 114

## Purpose

4Q ultrafast MaxQ NMR spectrum for a coupled four-spin system in the presence of realistic diffusion. Calculation time: hours, much faster on GPU

## Physical / mathematical content

- SPEN / ultrafast NMR examples. These files encode spatially dependent phase evolution and acquisition, linking pulse gradients, diffusion attenuation, and single-scan multidimensional encoding.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 12-13: Magnetic field; implemented by `sys.magnet=14.1`.
- Lines 15-16: Chemical shifts; implemented by `sys.isotopes={'1H','1H','1H','1H'}`.
- Lines 19-20: 3J couplings; implemented by `inter.coupling.scalar{1,2}=8.0`.
- Lines 24-25: 4J couplings; implemented by `inter.coupling.scalar{1,3}=3.0`.
- Lines 28-29: 5J couplings; implemented by `inter.coupling.scalar{1,4}=2.0`.
- Lines 32-33: Coherence selection; implemented by `parameters.mqorder=+4`.
- Lines 35-36: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 39-40: Algorithmic options; implemented by `sys.disable={'pt'}`.
- Lines 43-44: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 47-48: Sample geometry; implemented by `parameters.dims=0.015`.
- Lines 52-53: Relaxation phantom; implemented by `parameters.rlx_ph={}`.
- Lines 56-57: Initial and detection state phantoms; implemented by `parameters.rho0_ph={ones(parameters.npts,1)}`.
- Lines 62-63: Diffusion and flow; implemented by `parameters.u=zeros(parameters.npts,1)`.
- Lines 66-67: Acquisition parameters; implemented by `parameters.spins={'1H'}`.
- Lines 73-74: Calculate maximal k-value; implemented by `k_max=parameters.npoints/parameters.dims`.
- Lines 76-77: Get acquisition gradient duration; implemented by `Ta=parameters.deltat*parameters.npoints`.
- Lines 79-80: Get acquisition gradient amplitude; implemented by `parameters.Ga=((2*pi)*k_max)/(spin(parameters.spins{1})*Ta)`.
- Lines 82-83: Timings; implemented by `parameters.delay=0.041`.

### Key state/data transformations

- Lines 13: computes `sys.magnet` using `sys.magnet=14.1`.
- Lines 16: computes `sys.isotopes` using `sys.isotopes={'1H','1H','1H','1H'}`.
- Lines 17: computes `inter.zeeman.scalar` using `inter.zeeman.scalar={0.50,0.35,0.15,0}`.
- Lines 20: computes `inter.coupling.scalar{1,2}` using `inter.coupling.scalar{1,2}=8.0`.
- Lines 21: computes `inter.coupling.scalar{2,3}` using `inter.coupling.scalar{2,3}=8.0`.
- Lines 22: computes `inter.coupling.scalar{3,4}` using `inter.coupling.scalar{3,4}=8.0`.
- Lines 25: computes `inter.coupling.scalar{1,3}` using `inter.coupling.scalar{1,3}=3.0`.
- Lines 26: computes `inter.coupling.scalar{2,4}` using `inter.coupling.scalar{2,4}=3.0`.
- Lines 29: computes `inter.coupling.scalar{1,4}` using `inter.coupling.scalar{1,4}=2.0`.
- Lines 30: computes `inter.coupling.scalar{4,4}` using `inter.coupling.scalar{4,4}=0`.
- Lines 33: computes `parameters.mqorder` using `parameters.mqorder=+4`.
- Lines 36: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 37: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 40: computes `sys.disable` using `sys.disable={'pt'}`.
- Lines 41: computes `sys.enable` using `sys.enable={'greedy'}`.
- Lines 44: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 48: computes `parameters.dims` using `parameters.dims=0.015`.
- Lines 49: computes `parameters.npts` using `parameters.npts=500`.

## Implementation structure

- 4Q ultrafast MaxQ NMR spectrum for a coupled four-spin
- system in the presence of realistic diffusion.
- Calculation time: hours, much faster on GPU
- Magnetic field
- Chemical shifts
- 3J couplings
- 4J couplings
- 5J couplings
- Coherence selection
- Basis set
- Algorithmic options
- Spinach housekeeping

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `state()`, `spin()`, `imaging()`, `kfigure()`, `subplot()`, `fftshift()`, `plot_uf()`.
