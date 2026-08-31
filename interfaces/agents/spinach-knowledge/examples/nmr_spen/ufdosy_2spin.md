# examples/nmr_spen/ufdosy_2spin.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_spen/ufdosy_2spin.m`
- Signature: `ufdosy_2spin()`
- Total lines: 102

## Purpose

Ultrafast DOSY for two coupled spins with additional complications like DD and CSA relaxation, and spatial flow. Calculation time: minutes on NVidia Tesla A100, much longer on CPU Ludmilla Guduff Jean-Nicolas Dumez Ilya Kuprov

## Physical / mathematical content

- SPEN / ultrafast NMR examples. These files encode spatially dependent phase evolution and acquisition, linking pulse gradients, diffusion attenuation, and single-scan multidimensional encoding.
- The relaxation model is Redfield-type perturbation theory: fluctuating interactions enter through correlation functions or spectral densities and generate a linear relaxation superoperator.
- Chemical-shift anisotropy is present: shielding is treated as a second-rank tensor whose orientation relative to the field or rotor axis modulates line shapes and transfer dynamics.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.
- Finite-difference discretisation appears in the implementation, so numerical accuracy depends on stencil order, boundary handling, and the balance between resolution and conditioning.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 12-13: Spin system; implemented by `sys.isotopes={'1H','1H'}`.
- Lines 15-16: Interactions; implemented by `sys.magnet=14.1`.
- Lines 26-27: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 30-31: Relaxation theory; implemented by `inter.relaxation={'redfield'}`.
- Lines 36-37: Algorithmic options; implemented by `sys.disable={'pt'}`.
- Lines 40-41: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 44-45: Assumptions; implemented by `spin_system=assume(spin_system,'nmr')`.
- Lines 47-48: Sample geometry; implemented by `parameters.dims=0.015`.
- Lines 52-53: Relaxation phantom; implemented by `parameters.rlx_ph={ones(parameters.npts,1)}`.
- Lines 56-57: Initial and detection state phantoms; implemented by `parameters.rho0_ph={ones(parameters.npts,1)}`.
- Lines 62-63: Diffusion and flow; implemented by `parameters.u=1e-4*ones(parameters.npts,1)`.
- Lines 66-67: Acquisition parameters; implemented by `parameters.spins={'1H'}`.
- Lines 75-76: Encoding parameters; implemented by `parameters.pulsenpoints=1000`.
- Lines 84-85: Simulation; implemented by `fid=imaging(spin_system,@spendosy,parameters)`.
- Lines 87-88: Processing; implemented by `squarespectrum=fftshift(fft(fid,[],1),1)`.
- Lines 91-92: Plotting; implemented by `npoints2=size(squarespectrum,2)`.

### Key state/data transformations

- Lines 13: computes `sys.isotopes` using `sys.isotopes={'1H','1H'}`.
- Lines 16: computes `sys.magnet` using `sys.magnet=14.1`.
- Lines 17: computes `inter.zeeman.scalar` using `inter.zeeman.scalar={6.5,7.5}`.
- Lines 18: computes `inter.coupling.scalar{1,2}` using `inter.coupling.scalar{1,2}=15`.
- Lines 19: computes `inter.coupling.scalar{2,2}` using `inter.coupling.scalar{2,2}=0`.
- Lines 20: computes `inter.zeeman.eigs{1}` using `inter.zeeman.eigs{1}=[-10 -10 20]`.
- Lines 21: computes `inter.zeeman.euler{1}` using `inter.zeeman.euler{1}=[0 0 0]`.
- Lines 22: computes `inter.zeeman.eigs{2}` using `inter.zeeman.eigs{2}=[-10 -10 20]`.
- Lines 23: computes `inter.zeeman.euler{2}` using `inter.zeeman.euler{2}=[0 pi/2 0]`.
- Lines 24: computes `inter.coordinates` using `inter.coordinates={[0 0 0]; [0 0 2.5]}`.
- Lines 27: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 28: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 31: computes `inter.relaxation` using `inter.relaxation={'redfield'}`.
- Lines 32: computes `inter.tau_c` using `inter.tau_c={1.0e-9}`.
- Lines 33: computes `inter.rlx_keep` using `inter.rlx_keep='secular'`.
- Lines 34: computes `inter.equilibrium` using `inter.equilibrium='zero'`.
- Lines 37: computes `sys.disable` using `sys.disable={'pt'}`.
- Lines 38: computes `sys.enable` using `sys.enable={'greedy'}`.

## Implementation structure

- Ultrafast DOSY for two coupled spins with additional complications
- like DD and CSA relaxation, and spatial flow.
- Calculation time: minutes on NVidia Tesla A100, much longer on CPU
- Ludmilla Guduff
- Jean-Nicolas Dumez
- Ilya Kuprov
- Spin system
- Interactions
- Basis set
- Relaxation theory
- Algorithmic options
- Spinach housekeeping

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `assume()`, `relaxation()`, `state()`, `imaging()`, `fftshift()`, `spin()`, `kfigure()`, `kxlabel()`, `kylabel()`.
