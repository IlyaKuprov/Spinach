# examples/nmr_solids/hmqc_mas_dq.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_solids/hmqc_mas_dq.m`
- Signature: `hmqc_mas_dq()`
- Total lines: 74

## Purpose

Powder magic angle spinning CN2D experiment (rotor-synchronized de- tection) on a 14N-1H spin pair using 1D Fokker-Planck equation and a spherical grid. The calculation accounts for the second-order qu- adrupolar shift and lineshape. Calculation time: hours on CPU, minutes with a Tesla V100 GPU.

## Physical / mathematical content

- Solid-state NMR examples. The key physics is anisotropic spin interactions under static or magic-angle-spinning conditions: chemical-shift anisotropy, dipolar coupling, quadrupolar coupling, cross-polarisation, and orientation averaging using Floquet, Fokker-Planck, or direct powder quadrature formalisms.
- The file uses a Fokker-Planck-style enlarged state space in which spatial or orientational coordinates are promoted to extra dimensions and coupled to spin dynamics through differential operators.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 14-15: System specification; implemented by `sys.magnet=19.96`.
- Lines 23-24: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 27-31: Use GPU if present sys.enable={'gpu'};; implemented by `spin_system=create(sys,inter)`.
- Lines 30-31: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 34-35: Experiment setup; implemented by `parameters.rate=125000`.
- Lines 51-52: Simulation; implemented by `fid=singlerot(spin_system,@cn2d_dq,parameters,'qnmr')`.
- Lines 54-55: Apodisation; implemented by `fid.cos=apodisation(spin_system,fid.cos,{{'sqcos'},{'sqcos'}})`.
- Lines 58-59: F2 Fourier transform; implemented by `f1_cos=fftshift(fft(fid.cos,parameters.zerofill(2),1),1)`.
- Lines 62-63: Form States signal; implemented by `f1_states=real(f1_cos)+1i*real(f1_sin)`.
- Lines 65-66: F1 Fourier transform; implemented by `spectrum=fftshift(fft(f1_states,parameters.zerofill(1),2),2)`.
- Lines 68-69: Plotting; implemented by `kfigure(); scale_figure([1.5 2.0])`.

### Key state/data transformations

- Lines 15: computes `sys.magnet` using `sys.magnet=19.96`.
- Lines 16: computes `sys.isotopes` using `sys.isotopes={'14N','1H'}`.
- Lines 17: computes `inter.coupling.matrix{1,1}` using `inter.coupling.matrix{1,1}=eeqq2nqi(1.18e6,0.53,1,[0 0 0])`.
- Lines 18: computes `inter.coupling.matrix{2,2}` using `inter.coupling.matrix{2,2}=[]`.
- Lines 19: computes `inter.zeeman.scalar` using `inter.zeeman.scalar={32.4 5}`.
- Lines 20: computes `inter.coordinates` using `inter.coordinates={[0.00 0.00 0.00]`.
- Lines 24: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 25: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 31: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 35: computes `parameters.rate` using `parameters.rate=125000`.
- Lines 36: computes `parameters.axis` using `parameters.axis=[1 1 1]`.
- Lines 37: computes `parameters.max_rank` using `parameters.max_rank=16`.
- Lines 38: computes `parameters.grid` using `parameters.grid='rep_2ang_200pts_sph'`.
- Lines 39: computes `parameters.sweep` using `parameters.sweep=[parameters.rate/4 20000]`.
- Lines 40: computes `parameters.npoints` using `parameters.npoints=[256 128]`.
- Lines 41: computes `parameters.zerofill` using `parameters.zerofill=[1024 512]`.
- Lines 42: computes `parameters.spins` using `parameters.spins={'14N','1H'}`.
- Lines 43: computes `parameters.rframes` using `parameters.rframes={{'14N',3}}`.

## Implementation structure

- Powder magic angle spinning CN2D experiment (rotor-synchronized de-
- tection) on a 14N-1H spin pair using 1D Fokker-Planck equation and
- a spherical grid. The calculation accounts for the second-order qu-
- adrupolar shift and lineshape.
- Calculation time: hours on CPU, minutes with a Tesla V100 GPU.
- System specification
- Basis set
- Use GPU if present
- sys.enable={'gpu'};
- Spinach housekeeping
- Experiment setup
- Simulation

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `eeqq2nqi()`, `create()`, `basis()`, `state()`, `singlerot()`, `apodisation()`, `fftshift()`, `kfigure()`, `scale_figure()`, `plot_2d()`.
