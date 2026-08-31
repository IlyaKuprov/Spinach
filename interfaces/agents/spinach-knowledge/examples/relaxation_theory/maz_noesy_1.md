# examples/relaxation_theory/maz_noesy_1.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/relaxation_theory/maz_noesy_1.m`
- Signature: `maz_noesy_1()`
- Total lines: 137

## Purpose

15N-labelled methylaziridine NOESY, including the effects of the scalar relaxation of the first kind, caused by the modulation of J-coupling by the nitrogen centre inversion process. The calcula- tion illustrates the effect described in: Calculation time: minutes

## Physical / mathematical content

- Relaxation-theory examples. The mathematical backbone is Bloch-Redfield-Wangsness or stochastic Liouville theory, spectral densities, cross-correlation terms, motional models, and extraction of longitudinal/transverse decay behaviour from superoperators.
- The relaxation model is Redfield-type perturbation theory: fluctuating interactions enter through correlation functions or spectral densities and generate a linear relaxation superoperator.
- Propagation is accelerated with a Krylov-subspace method, replacing direct matrix exponentiation by projection into a much smaller Arnoldi/Lanczos-type subspace.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.
- A Krylov-subspace or Arnoldi construction is used to avoid forming or exponentiating very large dense propagators directly.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 16-17: Magnet induction; implemented by `sys.magnet=11.75`.
- Lines 19-20: Isotopes; implemented by `sys.isotopes={'1H','1H','1H','15N','1H','1H','1H','1H'}`.
- Lines 22-23: Absolute shielding (vacuum DFT); implemented by `inter.zeeman.matrix{1}=[32.7075 6.1977 -3.6061`.
- Lines 48-49: Assign isotropic components from the experiment; implemented by `expt_shift=[1.3 1.7 1.9 0.0 0.1 1.2 1.2 1.2]`.
- Lines 55-56: Scalar couplings (vacuum DFT); implemented by `inter.coupling.scalar=cell(8)`.
- Lines 68-69: Coordinates (Angstrom, vacuum DFT); implemented by `inter.coordinates={[-1.812977 -1.098554 0.444452]`.
- Lines 78-79: Algorithmic options; implemented by `sys.tols.inter_cutoff=2.0`.
- Lines 83-84: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 89-90: Relaxation superoperator; implemented by `inter.relaxation={'redfield','SRFK'}`.
- Lines 100-101: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 104-105: Pulse sequence parameters; implemented by `parameters.spins={'1H'}`.
- Lines 114-115: Simulation; implemented by `fid=liquid(spin_system,@noesy,parameters,'nmr')`.
- Lines 117-118: Apodisation; implemented by `fid.cos=apodisation(spin_system,fid.cos,{{'cos'},{'cos'}})`.
- Lines 121-122: F2 Fourier transform; implemented by `f1_cos=real(fftshift(fft(fid.cos,parameters.zerofill(2),1),1))`.
- Lines 125-126: States signal; implemented by `f1_states=f1_cos-1i*f1_sin`.
- Lines 128-129: F1 Fourier transform; implemented by `spectrum=fftshift(fft(f1_states,parameters.zerofill(1),2),2)`.
- Lines 131-132: Plotting; implemented by `kfigure(); scale_figure([1.5 2.0])`.

### Control flow inferred from the code

- Line 50: `for` loop over `n=1:numel(inter.zeeman.matrix)`.

### Key state/data transformations

- Lines 17: computes `sys.magnet` using `sys.magnet=11.75`.
- Lines 20: computes `sys.isotopes` using `sys.isotopes={'1H','1H','1H','15N','1H','1H','1H','1H'}`.
- Lines 23: computes `inter.zeeman.matrix{1}` using `inter.zeeman.matrix{1}=[32.7075 6.1977 -3.6061`.
- Lines 26: computes `inter.zeeman.matrix{2}` using `inter.zeeman.matrix{2}=[27.5505 -0.4604 -3.9670`.
- Lines 29: computes `inter.zeeman.matrix{3}` using `inter.zeeman.matrix{3}=[27.8325 -1.3662 -3.4664`.
- Lines 32: computes `inter.zeeman.matrix{4}` using `inter.zeeman.matrix{4}=[257.9722 16.9286 -28.3305`.
- Lines 35: computes `inter.zeeman.matrix{5}` using `inter.zeeman.matrix{5}=[28.3494 3.7981 -3.9745`.
- Lines 38: computes `inter.zeeman.matrix{6}` using `inter.zeeman.matrix{6}=[33.2953 4.3402 0.8099`.
- Lines 41: computes `inter.zeeman.matrix{7}` using `inter.zeeman.matrix{7}=[33.2302 -4.7306 1.2366`.
- Lines 44: computes `inter.zeeman.matrix{8}` using `inter.zeeman.matrix{8}=[31.4875 0.1078 -2.9339`.
- Lines 49: computes `expt_shift` using `expt_shift=[1.3 1.7 1.9 0.0 0.1 1.2 1.2 1.2]`.
- Lines 51: computes `inter.zeeman.matrix{n}` using `inter.zeeman.matrix{n}=inter.zeeman.matrix{n}-eye(3)*trace(inter.zeeman.matrix{n})/3`.
- Lines 56: computes `inter.coupling.scalar` using `inter.coupling.scalar=cell(8)`.
- Lines 57: computes `inter.coupling.scalar{1,2}` using `inter.coupling.scalar{1,2}=-11.6`.
- Lines 58: computes `inter.coupling.scalar{1,3}` using `inter.coupling.scalar{1,3}=5.5`.
- Lines 59: computes `inter.coupling.scalar{1,5}` using `inter.coupling.scalar{1,5}=7.5`.
- Lines 60: computes `inter.coupling.scalar{2,5}` using `inter.coupling.scalar{2,5}=12.7`.
- Lines 61: computes `inter.coupling.scalar{3,6}` using `inter.coupling.scalar{3,6}=4.5`.

## Implementation structure

- 15N-labelled methylaziridine NOESY, including the effects of the
- scalar relaxation of the first kind, caused by the modulation of
- J-coupling by the nitrogen centre inversion process. The calcula-
- tion illustrates the effect described in:
- Calculation time: minutes
- Magnet induction
- Isotopes
- Absolute shielding (vacuum DFT)
- Assign isotropic components from the experiment
- Scalar couplings (vacuum DFT)
- Coordinates (Angstrom, vacuum DFT)
- Algorithmic options

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `expt_shift()`, `create()`, `basis()`, `state()`, `liquid()`, `apodisation()`, `fftshift()`, `kfigure()`, `scale_figure()`, `plot_2d()`.
