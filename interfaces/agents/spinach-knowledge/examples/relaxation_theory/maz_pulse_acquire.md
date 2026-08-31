# examples/relaxation_theory/maz_pulse_acquire.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/relaxation_theory/maz_pulse_acquire.m`
- Signature: `maz_pulse_acquire()`
- Total lines: 132

## Purpose

Methylaziridine pulse-acquire, showing the effect of the scalar relaxation of the second kind due to the fast quadrupolar rela- xation of the 14N nuclei. Calculation time: minutes

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

- Lines 13-14: Magnet induction; implemented by `sys.magnet=11.75`.
- Lines 16-17: Isotopes; implemented by `sys.isotopes={'1H','1H','1H','14N','1H','1H','1H','1H'}`.
- Lines 19-20: Absolute shielding (vacuum DFT); implemented by `inter.zeeman.matrix{1}=[32.7075 6.1977 -3.6061`.
- Lines 45-46: Assign isotropic components from the experiment; implemented by `expt_shift=[1.3 1.7 1.9 0.0 0.1 1.2 1.2 1.2]`.
- Lines 52-53: Quadrupole couplings (vacuum DFT); implemented by `inter.coupling.matrix=cell(8)`.
- Lines 58-59: Scalar couplings (vacuum DFT); implemented by `inter.coupling.scalar=cell(8)`.
- Lines 75-76: Coordinates (Angstrom, vacuum DFT); implemented by `inter.coordinates={[-1.812977 -1.098554 0.444452]`.
- Lines 85-86: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 91-92: Relaxation superoperator; implemented by `inter.relaxation={'redfield','SRSK'}`.
- Lines 98-99: Algorithmic options; implemented by `sys.tols.inter_cutoff=2.0`.
- Lines 103-104: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 107-108: Sequence parameters; implemented by `parameters.spins={'1H'}`.
- Lines 119-120: Simulation; implemented by `fid=liquid(spin_system,@acquire,parameters,'nmr')`.
- Lines 122-123: Apodisation; implemented by `fid=apodisation(spin_system,fid,{{'exp',6}})`.
- Lines 125-126: Fourier transform; implemented by `spectrum=fftshift(fft(fid,parameters.zerofill))`.
- Lines 128-129: Plotting; implemented by `kfigure(); plot_1d(spin_system,real(spectrum),parameters)`.

### Control flow inferred from the code

- Line 47: `for` loop over `n=1:numel(inter.zeeman.matrix)`.

### Key state/data transformations

- Lines 14: computes `sys.magnet` using `sys.magnet=11.75`.
- Lines 17: computes `sys.isotopes` using `sys.isotopes={'1H','1H','1H','14N','1H','1H','1H','1H'}`.
- Lines 20: computes `inter.zeeman.matrix{1}` using `inter.zeeman.matrix{1}=[32.7075 6.1977 -3.6061`.
- Lines 23: computes `inter.zeeman.matrix{2}` using `inter.zeeman.matrix{2}=[27.5505 -0.4604 -3.9670`.
- Lines 26: computes `inter.zeeman.matrix{3}` using `inter.zeeman.matrix{3}=[27.8325 -1.3662 -3.4664`.
- Lines 29: computes `inter.zeeman.matrix{4}` using `inter.zeeman.matrix{4}=[257.9722 16.9286 -28.3305`.
- Lines 32: computes `inter.zeeman.matrix{5}` using `inter.zeeman.matrix{5}=[28.3494 3.7981 -3.9745`.
- Lines 35: computes `inter.zeeman.matrix{6}` using `inter.zeeman.matrix{6}=[33.2953 4.3402 0.8099`.
- Lines 38: computes `inter.zeeman.matrix{7}` using `inter.zeeman.matrix{7}=[33.2302 -4.7306 1.2366`.
- Lines 41: computes `inter.zeeman.matrix{8}` using `inter.zeeman.matrix{8}=[31.4875 0.1078 -2.9339`.
- Lines 46: computes `expt_shift` using `expt_shift=[1.3 1.7 1.9 0.0 0.1 1.2 1.2 1.2]`.
- Lines 48: computes `inter.zeeman.matrix{n}` using `inter.zeeman.matrix{n}=inter.zeeman.matrix{n}-eye(3)*trace(inter.zeeman.matrix{n})/3`.
- Lines 53: computes `inter.coupling.matrix` using `inter.coupling.matrix=cell(8)`.
- Lines 54: computes `inter.coupling.matrix{4,4}` using `inter.coupling.matrix{4,4}=[-1.2932 0.6251 1.8700`.
- Lines 59: computes `inter.coupling.scalar` using `inter.coupling.scalar=cell(8)`.
- Lines 60: computes `inter.coupling.scalar{1,2}` using `inter.coupling.scalar{1,2}=-11.6`.
- Lines 61: computes `inter.coupling.scalar{1,3}` using `inter.coupling.scalar{1,3}=5.5`.
- Lines 62: computes `inter.coupling.scalar{1,5}` using `inter.coupling.scalar{1,5}=7.5`.

## Implementation structure

- Methylaziridine pulse-acquire, showing the effect of the scalar
- relaxation of the second kind due to the fast quadrupolar rela-
- xation of the 14N nuclei.
- Calculation time: minutes
- Magnet induction
- Isotopes
- Absolute shielding (vacuum DFT)
- Assign isotropic components from the experiment
- Quadrupole couplings (vacuum DFT)
- Scalar couplings (vacuum DFT)
- Coordinates (Angstrom, vacuum DFT)
- Basis set

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `expt_shift()`, `create()`, `basis()`, `state()`, `liquid()`, `apodisation()`, `fftshift()`, `kfigure()`, `plot_1d()`.
