# examples/imaging/dpfgse_signal_select.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/imaging/dpfgse_signal_select.m`
- Signature: `dpfgse_signal_select()`
- Total lines: 99

## Purpose

DPFGSE signal selection example for a solution of GABA in water. Gradients and soft pulses are done explicitly. Simulation time: minutes, faster with a Tesla V100 GPU.

## Physical / mathematical content

- MRI and spectroscopic-imaging examples. These files combine gradient terms, spatial encoding, diffusion, slice selection, k-space sampling, and Fourier reconstruction, generally within Fokker-Planck or explicit spatial-grid descriptions.
- Propagation is accelerated with a Krylov-subspace method, replacing direct matrix exponentiation by projection into a much smaller Arnoldi/Lanczos-type subspace.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.
- A Krylov-subspace or Arnoldi construction is used to avoid forming or exponentiating very large dense propagators directly.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 12-13: Isotopes; implemented by `sys.isotopes={'1H','1H','1H','1H','1H','1H','1H'}`.
- Lines 15-16: Magnet (Tesla); implemented by `sys.magnet=5.9`.
- Lines 18-19: Chemical shifts (ppm); implemented by `inter.zeeman.scalar={3.00 3.00 1.88 1.88 2.28 2.28 4.80}`.
- Lines 21-22: J-couplings (Hz); implemented by `inter.coupling.scalar{1,3}=7.36`.
- Lines 32-33: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 38-39: Disable path tracing; implemented by `sys.disable={'pt','krylov'}`.
- Lines 41-45: This needs a GPU sys.enable={'gpu'};; implemented by `spin_system=create(sys,inter)`.
- Lines 44-45: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 48-49: Sequence parameters; implemented by `parameters.spins={'1H'}`.
- Lines 60-61: Signal selection; implemented by `pulse_nsteps=10`.
- Lines 68-69: Sample geometry; implemented by `parameters.dims=0.30`.
- Lines 73-74: Relaxation phantoms and operators; implemented by `parameters.rlx_ph={}; parameters.rlx_op={}`.
- Lines 76-77: Initial and detection state phantoms; implemented by `parameters.rho0_ph={ones(parameters.npts,1)}`.
- Lines 82-83: No diffusion or flow; implemented by `parameters.u=zeros(parameters.npts,1)`.
- Lines 86-87: Run simulation; implemented by `fid=imaging(spin_system,@dpfgse_select,parameters)`.
- Lines 89-90: Apodisation; implemented by `fid=apodisation(spin_system,fid,{{'exp',6}})`.
- Lines 92-93: Run the Fourier transform; implemented by `spectrum=fftshift(fft(fid,parameters.zerofill))`.
- Lines 95-96: Plotting; implemented by `kfigure(); plot_1d(spin_system,real(spectrum),parameters)`.

### Key state/data transformations

- Lines 13: computes `sys.isotopes` using `sys.isotopes={'1H','1H','1H','1H','1H','1H','1H'}`.
- Lines 16: computes `sys.magnet` using `sys.magnet=5.9`.
- Lines 19: computes `inter.zeeman.scalar` using `inter.zeeman.scalar={3.00 3.00 1.88 1.88 2.28 2.28 4.80}`.
- Lines 22: computes `inter.coupling.scalar{1,3}` using `inter.coupling.scalar{1,3}=7.36`.
- Lines 23: computes `inter.coupling.scalar{1,4}` using `inter.coupling.scalar{1,4}=7.36`.
- Lines 24: computes `inter.coupling.scalar{2,3}` using `inter.coupling.scalar{2,3}=7.36`.
- Lines 25: computes `inter.coupling.scalar{2,4}` using `inter.coupling.scalar{2,4}=7.36`.
- Lines 26: computes `inter.coupling.scalar{3,5}` using `inter.coupling.scalar{3,5}=7.58`.
- Lines 27: computes `inter.coupling.scalar{3,6}` using `inter.coupling.scalar{3,6}=7.58`.
- Lines 28: computes `inter.coupling.scalar{4,5}` using `inter.coupling.scalar{4,5}=7.58`.
- Lines 29: computes `inter.coupling.scalar{4,6}` using `inter.coupling.scalar{4,6}=7.58`.
- Lines 30: computes `inter.coupling.scalar{7,7}` using `inter.coupling.scalar{7,7}=0`.
- Lines 33: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 34: computes `bas.approximation` using `bas.approximation='IK-2'`.
- Lines 35: computes `bas.space_level` using `bas.space_level=1`.
- Lines 36: computes `bas.connectivity` using `bas.connectivity='scalar_couplings'`.
- Lines 39: computes `sys.disable` using `sys.disable={'pt','krylov'}`.
- Lines 45: computes `spin_system` using `spin_system=create(sys,inter)`.

## Implementation structure

- DPFGSE signal selection example for a solution of GABA
- in water. Gradients and soft pulses are done explicitly.
- Simulation time: minutes, faster with a Tesla V100 GPU.
- Isotopes
- Magnet (Tesla)
- Chemical shifts (ppm)
- J-couplings (Hz)
- Basis set
- Disable path tracing
- This needs a GPU
- sys.enable={'gpu'};
- Spinach housekeeping

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `pulse_shape()`, `state()`, `imaging()`, `apodisation()`, `fftshift()`, `kfigure()`, `plot_1d()`.
