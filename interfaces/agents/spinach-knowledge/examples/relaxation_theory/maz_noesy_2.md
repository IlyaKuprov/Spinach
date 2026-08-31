# examples/relaxation_theory/maz_noesy_2.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/relaxation_theory/maz_noesy_2.m`
- Signature: `maz_noesy_2()`
- Total lines: 145

## Purpose

Methylaziridine NOESY, including the effects of the scalar relaxation of the first kind (caused by the modulation of J-coupling by the nit- rogen centre inversion process) and second kind (caused by the rapid quadrupolar relaxation of the 14N nucleus). The calculation illustra- tes the effect described in: Calculation time: minutes.

## Physical / mathematical content

- Relaxation-theory examples. The mathematical backbone is Bloch-Redfield-Wangsness or stochastic Liouville theory, spectral densities, cross-correlation terms, motional models, and extraction of longitudinal/transverse decay behaviour from superoperators.
- The relaxation model is Redfield-type perturbation theory: fluctuating interactions enter through correlation functions or spectral densities and generate a linear relaxation superoperator.
- Propagation is accelerated with a Krylov-subspace method, replacing direct matrix exponentiation by projection into a much smaller Arnoldi/Lanczos-type subspace.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.
- Quadrupolar physics is relevant: nuclei with spin > 1/2 interact with the electric field gradient tensor, introducing second-rank anisotropy, asymmetry, and overtone or MQ phenomena.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.
- A Krylov-subspace or Arnoldi construction is used to avoid forming or exponentiating very large dense propagators directly.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 17-18: Magnet induction; implemented by `sys.magnet=11.75`.
- Lines 20-21: Isotopes; implemented by `sys.isotopes={'1H','1H','1H','14N','1H','1H','1H','1H'}`.
- Lines 23-24: Absolute shielding tensors (vacuum DFT); implemented by `inter.zeeman.matrix{1}=[32.7075 6.1977 -3.6061`.
- Lines 49-50: Assign isotropic components from the experiment; implemented by `expt_shift=[1.3 1.7 1.9 0.0 0.1 1.2 1.2 1.2]`.
- Lines 56-57: Quadrupole couplings (vacuum DFT); implemented by `inter.coupling.matrix=cell(8)`.
- Lines 62-63: Scalar couplings (vacuum DFT); implemented by `inter.coupling.scalar=cell(8)`.
- Lines 75-76: Coordinates (Angstrom, vacuum DFT); implemented by `inter.coordinates={[-1.812977 -1.098554 0.444452]`.
- Lines 85-86: Algorithmic options; implemented by `sys.tols.inter_cutoff=2.0`.
- Lines 90-91: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 96-97: Relaxation superoperator; implemented by `inter.relaxation={'redfield','SRFK','SRSK'}`.
- Lines 108-109: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 112-113: Sequence parameters; implemented by `parameters.spins={'1H'}`.
- Lines 122-123: Simulation; implemented by `fid=liquid(spin_system,@noesy,parameters,'nmr')`.
- Lines 125-126: Apodisation; implemented by `fid.cos=apodisation(spin_system,fid.cos,{{'cos'},{'cos'}})`.
- Lines 129-130: F2 Fourier transform; implemented by `f1_cos=real(fftshift(fft(fid.cos,parameters.zerofill(2),1),1))`.
- Lines 133-134: States signal; implemented by `f1_states=f1_cos-1i*f1_sin`.
- Lines 136-137: F1 Fourier transform; implemented by `spectrum=fftshift(fft(f1_states,parameters.zerofill(1),2),2)`.
- Lines 139-140: Plotting; implemented by `kfigure(); scale_figure([1.5 2.0])`.

### Control flow inferred from the code

- Line 51: `for` loop over `n=1:numel(inter.zeeman.matrix)`.

### Key state/data transformations

- Lines 18: computes `sys.magnet` using `sys.magnet=11.75`.
- Lines 21: computes `sys.isotopes` using `sys.isotopes={'1H','1H','1H','14N','1H','1H','1H','1H'}`.
- Lines 24: computes `inter.zeeman.matrix{1}` using `inter.zeeman.matrix{1}=[32.7075 6.1977 -3.6061`.
- Lines 27: computes `inter.zeeman.matrix{2}` using `inter.zeeman.matrix{2}=[27.5505 -0.4604 -3.9670`.
- Lines 30: computes `inter.zeeman.matrix{3}` using `inter.zeeman.matrix{3}=[27.8325 -1.3662 -3.4664`.
- Lines 33: computes `inter.zeeman.matrix{4}` using `inter.zeeman.matrix{4}=[257.9722 16.9286 -28.3305`.
- Lines 36: computes `inter.zeeman.matrix{5}` using `inter.zeeman.matrix{5}=[28.3494 3.7981 -3.9745`.
- Lines 39: computes `inter.zeeman.matrix{6}` using `inter.zeeman.matrix{6}=[33.2953 4.3402 0.8099`.
- Lines 42: computes `inter.zeeman.matrix{7}` using `inter.zeeman.matrix{7}=[33.2302 -4.7306 1.2366`.
- Lines 45: computes `inter.zeeman.matrix{8}` using `inter.zeeman.matrix{8}=[31.4875 0.1078 -2.9339`.
- Lines 50: computes `expt_shift` using `expt_shift=[1.3 1.7 1.9 0.0 0.1 1.2 1.2 1.2]`.
- Lines 52: computes `inter.zeeman.matrix{n}` using `inter.zeeman.matrix{n}=inter.zeeman.matrix{n}-eye(3)*trace(inter.zeeman.matrix{n})/3`.
- Lines 57: computes `inter.coupling.matrix` using `inter.coupling.matrix=cell(8)`.
- Lines 58: computes `inter.coupling.matrix{4,4}` using `inter.coupling.matrix{4,4}=[-1.2932 0.6251 1.8700`.
- Lines 63: computes `inter.coupling.scalar` using `inter.coupling.scalar=cell(8)`.
- Lines 64: computes `inter.coupling.scalar{1,2}` using `inter.coupling.scalar{1,2}=-11.6`.
- Lines 65: computes `inter.coupling.scalar{1,3}` using `inter.coupling.scalar{1,3}=5.5`.
- Lines 66: computes `inter.coupling.scalar{1,5}` using `inter.coupling.scalar{1,5}=7.5`.

## Implementation structure

- Methylaziridine NOESY, including the effects of the scalar relaxation
- of the first kind (caused by the modulation of J-coupling by the nit-
- rogen centre inversion process) and second kind (caused by the rapid
- quadrupolar relaxation of the 14N nucleus). The calculation illustra-
- tes the effect described in:
- Calculation time: minutes.
- Magnet induction
- Isotopes
- Absolute shielding tensors (vacuum DFT)
- Assign isotropic components from the experiment
- Quadrupole couplings (vacuum DFT)
- Scalar couplings (vacuum DFT)

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `expt_shift()`, `create()`, `basis()`, `state()`, `liquid()`, `apodisation()`, `fftshift()`, `kfigure()`, `scale_figure()`, `plot_2d()`.
