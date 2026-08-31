# examples/nmr_spen/ufmq_6spin.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_spen/ufmq_6spin.m`
- Signature: `ufmq_6spin()`
- Total lines: 125

## Purpose

6Q ultrafast MaxQ NMR spectrum for a coupled six-spin system in the presence of realistic diffusion. Calculation time: hours on NVidia Tesla A100, much longer on CPU

## Physical / mathematical content

- SPEN / ultrafast NMR examples. These files encode spatially dependent phase evolution and acquisition, linking pulse gradients, diffusion attenuation, and single-scan multidimensional encoding.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 12-13: Magnetic field; implemented by `sys.magnet=14.1`.
- Lines 15-16: Chemical shifts; implemented by `sys.isotopes={'1H','1H','1H','1H','1H','1H'}`.
- Lines 19-20: 3J couplings; implemented by `inter.coupling.scalar{1,2}=8.00`.
- Lines 25-26: 4J couplings; implemented by `inter.coupling.scalar{1,3}=4.00`.
- Lines 33-34: 5J couplings; implemented by `inter.coupling.scalar{2,5}=4.00`.
- Lines 38-39: 6J couplings; implemented by `inter.coupling.scalar{1,5}=2.00`.
- Lines 43-44: Coherence selection; implemented by `parameters.mqorder=+6`.
- Lines 46-47: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 50-51: Algorithmic options; implemented by `sys.disable={'pt'}`.
- Lines 54-55: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 58-59: Sample geometry; implemented by `parameters.dims=0.015`.
- Lines 63-64: Relaxation phantom; implemented by `parameters.rlx_ph={}`.
- Lines 67-68: Initial and detection state phantoms; implemented by `parameters.rho0_ph={ones(parameters.npts,1)}`.
- Lines 73-74: Diffusion and flow; implemented by `parameters.u=zeros(parameters.npts,1)`.
- Lines 77-78: Acquisition parameters; implemented by `parameters.spins={'1H'}`.
- Lines 84-85: Calculate maximal k-value; implemented by `k_max=parameters.npoints/parameters.dims`.
- Lines 87-88: Get acquisition gradient duration; implemented by `Ta=parameters.deltat*parameters.npoints`.
- Lines 90-91: Get acquisition gradient amplitude; implemented by `parameters.Ga=((2*pi)*k_max)/(spin(parameters.spins{1})*Ta)`.

### Key state/data transformations

- Lines 13: computes `sys.magnet` using `sys.magnet=14.1`.
- Lines 16: computes `sys.isotopes` using `sys.isotopes={'1H','1H','1H','1H','1H','1H'}`.
- Lines 17: computes `inter.zeeman.scalar` using `inter.zeeman.scalar={-1.0 -0.5 0.0 +0.3 +0.7 +0.9}`.
- Lines 20: computes `inter.coupling.scalar{1,2}` using `inter.coupling.scalar{1,2}=8.00`.
- Lines 21: computes `inter.coupling.scalar{2,3}` using `inter.coupling.scalar{2,3}=8.00`.
- Lines 22: computes `inter.coupling.scalar{3,4}` using `inter.coupling.scalar{3,4}=8.00`.
- Lines 23: computes `inter.coupling.scalar{4,5}` using `inter.coupling.scalar{4,5}=8.00`.
- Lines 26: computes `inter.coupling.scalar{1,3}` using `inter.coupling.scalar{1,3}=4.00`.
- Lines 27: computes `inter.coupling.scalar{1,4}` using `inter.coupling.scalar{1,4}=4.00`.
- Lines 28: computes `inter.coupling.scalar{2,4}` using `inter.coupling.scalar{2,4}=4.00`.
- Lines 29: computes `inter.coupling.scalar{3,5}` using `inter.coupling.scalar{3,5}=4.00`.
- Lines 30: computes `inter.coupling.scalar{1,6}` using `inter.coupling.scalar{1,6}=4.00`.
- Lines 31: computes `inter.coupling.scalar{5,6}` using `inter.coupling.scalar{5,6}=4.00`.
- Lines 34: computes `inter.coupling.scalar{2,5}` using `inter.coupling.scalar{2,5}=4.00`.
- Lines 35: computes `inter.coupling.scalar{2,6}` using `inter.coupling.scalar{2,6}=4.00`.
- Lines 36: computes `inter.coupling.scalar{4,6}` using `inter.coupling.scalar{4,6}=4.00`.
- Lines 39: computes `inter.coupling.scalar{1,5}` using `inter.coupling.scalar{1,5}=2.00`.
- Lines 40: computes `inter.coupling.scalar{3,6}` using `inter.coupling.scalar{3,6}=2.00`.

## Implementation structure

- 6Q ultrafast MaxQ NMR spectrum for a coupled six-spin
- system in the presence of realistic diffusion.
- Calculation time: hours on NVidia Tesla A100, much longer on CPU
- Magnetic field
- Chemical shifts
- 3J couplings
- 4J couplings
- 5J couplings
- 6J couplings
- Coherence selection
- Basis set
- Algorithmic options

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `state()`, `spin()`, `imaging()`, `figure()`, `subplot()`, `fftshift()`, `plot_uf()`.
