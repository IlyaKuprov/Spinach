# examples/kinetics/aziridine_exsy_2.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/kinetics/aziridine_exsy_2.m`
- Signature: `aziridine_exsy_2()`
- Total lines: 205

## Purpose

NOESY/EXSY experiment on phenylaziridine, including scalar relaxation of the second kind induced by the 14N nucleus, in a situation where the chemical exchange is intermediate, lines are broadened and scalar relaxation of the first kind must be accounted for. Set to reproduce Figures 1a and 4a from All parameters, except for the isotropic chemical shifts, exchange ra- tes and correlation times come from a DFT calcula

## Physical / mathematical content

- Chemical-kinetics examples. The files couple spin dynamics to exchange, pumping, or nonlinear reaction networks represented by kinetic generators in Liouville space.
- The relaxation model is Redfield-type perturbation theory: fluctuating interactions enter through correlation functions or spectral densities and generate a linear relaxation superoperator.
- Propagation is accelerated with a Krylov-subspace method, replacing direct matrix exponentiation by projection into a much smaller Arnoldi/Lanczos-type subspace.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.
- Quadrupolar physics is relevant: nuclei with spin > 1/2 interact with the electric field gradient tensor, introducing second-rank anisotropy, asymmetry, and overtone or MQ phenomena.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.
- A Krylov-subspace or Arnoldi construction is used to avoid forming or exponentiating very large dense propagators directly.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 20-21: Magnet induction; implemented by `sys.magnet=11.75`.
- Lines 23-25: Isotopes; implemented by `sys.isotopes={'1H','1H','1H','14N','1H','1H','1H','1H','1H','1H', '1H','1H','1H','14N','1H','1H','1H','1H','1H','1H'}`.
- Lines 27-28: Coordinates (Angstrom); implemented by `inter.coordinates={[-3.41285300 0.73156000 1.07036600]`.
- Lines 49-50: 14N quadrupolar coupling; implemented by `inter.coupling.matrix=cell(20)`.
- Lines 58-59: Scalar couplings; implemented by `inter.coupling.scalar=cell(20)`.
- Lines 99-100: Chemical shift tensors; implemented by `inter.zeeman.matrix=cell(1,20)`.
- Lines 122-124: Experimental chemical shifts; implemented by `inter.zeeman.matrix=shift_iso(inter.zeeman.matrix,1:20,[2.237 1.576 2.970 0.0 1.134 7.250 7.250 7.250 7.320 7.320 2.060 2.111 2.843 0.0 0.536 7.250 7.250 7.250 7.320 7.3…`.
- Lines 125-126: Chemical kinetics; implemented by `kplus=1.2e3; kminus=1.2e3`.
- Lines 133-134: Relaxation theory; implemented by `inter.relaxation={'redfield','SRFK','SRSK'}`.
- Lines 152-153: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 159-160: Disable Krylov algorithm; implemented by `sys.disable={'krylov'}`.
- Lines 163-164: Proximity cut-off; implemented by `sys.tols.prox_cutoff=10.0`.
- Lines 166-167: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 170-171: Sequence parameters; implemented by `parameters.offset=2000`.
- Lines 179-180: Concentration-aware initial state; implemented by `parameters.rho0=state(spin_system,'Lz','1H','chem')`.
- Lines 182-183: Simulation; implemented by `fid=liquid(spin_system,@noesy,parameters,'nmr')`.
- Lines 185-186: Apodisation; implemented by `fid.cos=apodisation(spin_system,fid.cos,{{'sqcos'},{'sqcos'}})`.
- Lines 189-190: F2 Fourier transform; implemented by `f1_cos=real(fftshift(fft(fid.cos,parameters.zerofill(2),1),1))`.

### Key state/data transformations

- Lines 21: computes `sys.magnet` using `sys.magnet=11.75`.
- Lines 24-25: computes `sys.isotopes` using `sys.isotopes={'1H','1H','1H','14N','1H','1H','1H','1H','1H','1H', '1H','1H','1H','14N','1H','1H','1H','1H','1H','1H'}`.
- Lines 28: computes `inter.coordinates` using `inter.coordinates={[-3.41285300 0.73156000 1.07036600]`.
- Lines 50: computes `inter.coupling.matrix` using `inter.coupling.matrix=cell(20)`.
- Lines 51: computes `inter.coupling.matrix{4,4}` using `inter.coupling.matrix{4,4}=1e6*[-2.108e+000 1.299e+000 -1.434e+000`.
- Lines 54: computes `inter.coupling.matrix{14,14}` using `inter.coupling.matrix{14,14}=1e6*[-2.108e+000 1.299e+000 -1.434e+000`.
- Lines 59: computes `inter.coupling.scalar` using `inter.coupling.scalar=cell(20)`.
- Lines 60: computes `inter.coupling.scalar{2,1}` using `inter.coupling.scalar{2,1}= -1.122e+001`.
- Lines 61: computes `inter.coupling.scalar{3,1}` using `inter.coupling.scalar{3,1}= 5.824e+000`.
- Lines 62: computes `inter.coupling.scalar{3,2}` using `inter.coupling.scalar{3,2}= 1.867e+000`.
- Lines 63: computes `inter.coupling.scalar{4,1}` using `inter.coupling.scalar{4,1}= 4.404e+000`.
- Lines 64: computes `inter.coupling.scalar{4,2}` using `inter.coupling.scalar{4,2}= -5.709e-001`.
- Lines 65: computes `inter.coupling.scalar{4,3}` using `inter.coupling.scalar{4,3}= 5.049e+000`.
- Lines 66: computes `inter.coupling.scalar{5,1}` using `inter.coupling.scalar{5,1}= 7.729e+000`.
- Lines 67: computes `inter.coupling.scalar{5,2}` using `inter.coupling.scalar{5,2}= 1.308e+001`.
- Lines 68: computes `inter.coupling.scalar{5,3}` using `inter.coupling.scalar{5,3}= 8.069e+000`.
- Lines 69: computes `inter.coupling.scalar{6,3}` using `inter.coupling.scalar{6,3}= -1.481e+000`.
- Lines 70: computes `inter.coupling.scalar{7,3}` using `inter.coupling.scalar{7,3}= -1.548e+000`.

## Implementation structure

- NOESY/EXSY experiment on phenylaziridine, including scalar relaxation
- of the second kind induced by the 14N nucleus, in a situation where
- the chemical exchange is intermediate, lines are broadened and scalar
- relaxation of the first kind must be accounted for. Set to reproduce
- Figures 1a and 4a from
- All parameters, except for the isotropic chemical shifts, exchange ra-
- tes and correlation times come from a DFT calculation.
- Calculation time: hours, faster on a GPU.
- Magnet induction
- Isotopes
- Coordinates (Angstrom)
- 14N quadrupolar coupling

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `shift_iso()`, `create()`, `basis()`, `state()`, `liquid()`, `apodisation()`, `fftshift()`, `kfigure()`, `scale_figure()`, `plot_2d()`.
