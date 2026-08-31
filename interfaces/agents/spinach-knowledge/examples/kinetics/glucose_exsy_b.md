# examples/kinetics/glucose_exsy_b.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/kinetics/glucose_exsy_b.m`
- Signature: `glucose_exsy_b()`
- Total lines: 137

## Purpose

2D EXSY of transmembrane exchange of 3,3-difluoroglucose. See the fitting example set for the script that yielded the para- meters used below. Calculation time: seconds

## Physical / mathematical content

- Chemical-kinetics examples. The files couple spin dynamics to exchange, pumping, or nonlinear reaction networks represented by kinetic generators in Liouville space.
- The relaxation model is Redfield-type perturbation theory: fluctuating interactions enter through correlation functions or spectral densities and generate a linear relaxation superoperator.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 13-14: Magnet field; implemented by `sys.magnet=9.3933`.
- Lines 16-20: Isotopes; implemented by `sys.isotopes={'19F','19F', '19F','19F', '19F','19F', '19F','19F'}`.
- Lines 22-23: Chemical shifts; implemented by `inter.zeeman.scalar={-113.6784 -129.6771`.
- Lines 28-29: J-couplings; implemented by `inter.coupling.scalar=cell(8,8)`.
- Lines 35-36: Cartesian coordinates; implemented by `inter.coordinates={[-0.0551 -1.2087 -1.6523]`.
- Lines 48-49: Chemical subsystems; implemented by `inter.chem.parts={[1 2],`.
- Lines 54-55: Reaction rate matrix; implemented by `inter.chem.rates=[-0.3438 0.1550 0 0`.
- Lines 60-61: Equilibrium concentrations with alpha-beta imbalance; implemented by `inter.chem.concs=equilibrate(inter.chem.rates,[3.8034; 0; 14.2442; 0])`.
- Lines 63-64: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 67-68: Relaxation theory parameters; implemented by `inter.relaxation={'redfield','t1_t2'}`.
- Lines 78-79: Do not draw colorbars; implemented by `sys.disable={'colorbar'}`.
- Lines 81-82: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 85-86: Sequence parameters; implemented by `parameters.tmix=0.5`.
- Lines 95-96: Simulation; implemented by `fid=liquid(spin_system,@noesy,parameters,'nmr')`.
- Lines 98-99: Apodisation; implemented by `fid.cos=apodisation(spin_system,fid.cos,{{'sqcos'},{'sqcos'}})`.
- Lines 102-103: F2 Fourier transform; implemented by `f1_cos=real(fftshift(fft(fid.cos,parameters.zerofill(2),1),1))`.
- Lines 106-107: States signal; implemented by `f1_states=f1_cos-1i*f1_sin`.
- Lines 109-110: F1 Fourier transform; implemented by `spectrum=fftshift(fft(f1_states,parameters.zerofill(1),2),2)`.

### Key state/data transformations

- Lines 14: computes `sys.magnet` using `sys.magnet=9.3933`.
- Lines 17-20: computes `sys.isotopes` using `sys.isotopes={'19F','19F', '19F','19F', '19F','19F', '19F','19F'}`.
- Lines 23: computes `inter.zeeman.scalar` using `inter.zeeman.scalar={-113.6784 -129.6771`.
- Lines 29: computes `inter.coupling.scalar` using `inter.coupling.scalar=cell(8,8)`.
- Lines 30: computes `inter.coupling.scalar{1,2}` using `inter.coupling.scalar{1,2}=238.0633`.
- Lines 31: computes `inter.coupling.scalar{3,4}` using `inter.coupling.scalar{3,4}=238.0633`.
- Lines 32: computes `inter.coupling.scalar{5,6}` using `inter.coupling.scalar{5,6}=239.2315`.
- Lines 33: computes `inter.coupling.scalar{7,8}` using `inter.coupling.scalar{7,8}=239.2315`.
- Lines 36: computes `inter.coordinates` using `inter.coordinates={[-0.0551 -1.2087 -1.6523]`.
- Lines 49: computes `inter.chem.parts` using `inter.chem.parts={[1 2],`.
- Lines 55: computes `inter.chem.rates` using `inter.chem.rates=[-0.3438 0.1550 0 0`.
- Lines 61: computes `inter.chem.concs` using `inter.chem.concs=equilibrate(inter.chem.rates,[3.8034; 0; 14.2442; 0])`.
- Lines 64: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 65: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 68: computes `inter.relaxation` using `inter.relaxation={'redfield','t1_t2'}`.
- Lines 69: computes `inter.equilibrium` using `inter.equilibrium='zero'`.
- Lines 70: computes `inter.rlx_keep` using `inter.rlx_keep='secular'`.
- Lines 71: computes `inter.tau_c` using `inter.tau_c={0.9601e-9`.

## Implementation structure

- 2D EXSY of transmembrane exchange of 3,3-difluoroglucose. See
- the fitting example set for the script that yielded the para-
- meters used below.
- Calculation time: seconds
- Magnet field
- Isotopes
- Chemical shifts
- J-couplings
- Cartesian coordinates
- Chemical subsystems
- Reaction rate matrix
- Equilibrium concentrations with alpha-beta imbalance

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `equilibrate()`, `num2cell()`, `create()`, `basis()`, `state()`, `liquid()`, `apodisation()`, `fftshift()`, `load()`, `keep_rank()`, `atranspose()`, `kfigure()`, `scale_figure()`, `subplot()`, `plot_2d()`, `ktitle()`.
