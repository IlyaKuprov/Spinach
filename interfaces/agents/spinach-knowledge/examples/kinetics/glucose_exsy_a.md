# examples/kinetics/glucose_exsy_a.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/kinetics/glucose_exsy_a.m`
- Signature: `glucose_exsy_a()`
- Total lines: 168

## Purpose

2D EXSY of transmembrane exchange of 2,2,3,3-tetrafluoroglucose. See the fitting example set for the script that yielded the parameters used below. Calculation time: seconds

## Physical / mathematical content

- Chemical-kinetics examples. The files couple spin dynamics to exchange, pumping, or nonlinear reaction networks represented by kinetic generators in Liouville space.
- The relaxation model is Redfield-type perturbation theory: fluctuating interactions enter through correlation functions or spectral densities and generate a linear relaxation superoperator.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 13-14: Magnet field; implemented by `sys.magnet=9.4`.
- Lines 16-20: Isotopes; implemented by `sys.isotopes={'19F','19F','19F','19F', '19F','19F','19F','19F', '19F','19F','19F','19F', '19F','19F','19F','19F'}`.
- Lines 22-23: Chemical shifts; implemented by `inter.zeeman.scalar={-120.5380 -133.9429 -129.3169 -129.5320`.
- Lines 28-29: J-couplings; implemented by `inter.coupling.scalar=cell(16,16)`.
- Lines 59-60: Cartesian coordinates; implemented by `inter.coordinates={[-0.0551 -1.2087 -1.6523]`.
- Lines 80-81: Chemical subsystems; implemented by `inter.chem.parts={[1 2 3 4],`.
- Lines 86-87: Reaction rate matrix; implemented by `inter.chem.rates=[-1.0045 1.7738 0 0`.
- Lines 92-93: Equilibrate translocation with alpha-beta imbalance as the start; implemented by `inter.chem.concs=equilibrate(inter.chem.rates,[3.2258; 0; 3.1902; 0])`.
- Lines 95-96: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 99-100: Relaxation theory parameters; implemented by `inter.relaxation={'redfield'}`.
- Lines 108-109: Do not draw colorbars; implemented by `sys.disable={'colorbar'}`.
- Lines 111-112: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 115-116: Sequence parameters; implemented by `parameters.tmix=0.5`.
- Lines 125-126: Simulation; implemented by `fid=liquid(spin_system,@noesy,parameters,'nmr')`.
- Lines 128-129: Apodisation; implemented by `fid.cos=apodisation(spin_system,fid.cos,{{'sqcos'},{'sqcos'}})`.
- Lines 132-133: F2 Fourier transform; implemented by `f1_cos=real(fftshift(fft(fid.cos,parameters.zerofill(2),1),1))`.
- Lines 136-137: States signal; implemented by `f1_states=f1_cos-1i*f1_sin`.
- Lines 139-140: F1 Fourier transform; implemented by `spectrum=fftshift(fft(f1_states,parameters.zerofill(1),2),2)`.

### Key state/data transformations

- Lines 14: computes `sys.magnet` using `sys.magnet=9.4`.
- Lines 17-20: computes `sys.isotopes` using `sys.isotopes={'19F','19F','19F','19F', '19F','19F','19F','19F', '19F','19F','19F','19F', '19F','19F','19F','19F'}`.
- Lines 23: computes `inter.zeeman.scalar` using `inter.zeeman.scalar={-120.5380 -133.9429 -129.3169 -129.5320`.
- Lines 29: computes `inter.coupling.scalar` using `inter.coupling.scalar=cell(16,16)`.
- Lines 31: computes `inter.coupling.scalar{1,2}` using `inter.coupling.scalar{1,2}= 271.2924`.
- Lines 32: computes `inter.coupling.scalar{3,4}` using `inter.coupling.scalar{3,4}= 271.2924`.
- Lines 33: computes `inter.coupling.scalar{1,3}` using `inter.coupling.scalar{1,3}= 0.5401`.
- Lines 34: computes `inter.coupling.scalar{1,4}` using `inter.coupling.scalar{1,4}= -25.9884`.
- Lines 35: computes `inter.coupling.scalar{2,3}` using `inter.coupling.scalar{2,3}= 9.9625`.
- Lines 36: computes `inter.coupling.scalar{2,4}` using `inter.coupling.scalar{2,4}= -40.7675`.
- Lines 38: computes `inter.coupling.scalar{5,6}` using `inter.coupling.scalar{5,6}= 271.2924`.
- Lines 39: computes `inter.coupling.scalar{7,8}` using `inter.coupling.scalar{7,8}= 271.2924`.
- Lines 40: computes `inter.coupling.scalar{5,7}` using `inter.coupling.scalar{5,7}= 0.5401`.
- Lines 41: computes `inter.coupling.scalar{5,8}` using `inter.coupling.scalar{5,8}= -25.9884`.
- Lines 42: computes `inter.coupling.scalar{6,7}` using `inter.coupling.scalar{6,7}= 9.9625`.
- Lines 43: computes `inter.coupling.scalar{6,8}` using `inter.coupling.scalar{6,8}= -40.7675`.
- Lines 45: computes `inter.coupling.scalar{9,10}` using `inter.coupling.scalar{9,10}= 263.2659`.
- Lines 46: computes `inter.coupling.scalar{11,12}` using `inter.coupling.scalar{11,12}= 263.2659`.

## Implementation structure

- 2D EXSY of transmembrane exchange of 2,2,3,3-tetrafluoroglucose. See
- the fitting example set for the script that yielded the parameters
- used below.
- Calculation time: seconds
- Magnet field
- Isotopes
- Chemical shifts
- J-couplings
- Cartesian coordinates
- Chemical subsystems
- Reaction rate matrix
- Equilibrate translocation with alpha-beta imbalance as the start

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `equilibrate()`, `create()`, `basis()`, `state()`, `liquid()`, `apodisation()`, `fftshift()`, `load()`, `rot90()`, `keep_rank()`, `kfigure()`, `scale_figure()`, `subplot()`, `plot_2d()`, `ktitle()`, `histogram()`.
