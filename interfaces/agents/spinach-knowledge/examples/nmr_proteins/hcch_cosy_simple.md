# examples/nmr_proteins/hcch_cosy_simple.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_proteins/hcch_cosy_simple.m`
- Signature: `hcch_cosy_simple()`
- Total lines: 81

## Purpose

3D HCCH COSY experiment on a small protein fragment. Calculation time: seconds.

## Physical / mathematical content

- Protein NMR examples. These files specialise liquid-state pulse sequences to labelled biomolecules, exploiting one-bond and two-bond heteronuclear couplings, coherence pathway filtering, selective decoupling, and high-dimensional indirect detection.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 10-11: Spin system; implemented by `sys.isotopes={'13C','13C','13C','1H','1H'}`.
- Lines 14-15: Interactions; implemented by `sys.magnet=14.1`.
- Lines 27-28: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 31-32: Sequence parameters; implemented by `parameters.J_ch=140`.
- Lines 43-44: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 47-48: Simulation; implemented by `fid=liquid(spin_system,@hcch_cosy,parameters,'nmr')`.
- Lines 50-51: Apodisation; implemented by `fid.pos_pos=apodisation(spin_system,fid.pos_pos,{{'sqcos'},{'sqcos'},{'sqcos'}})`.
- Lines 56-57: F3 Fourier transform; implemented by `f3_pos_pos=fftshift(fft(fid.pos_pos,parameters.zerofill(3),3),3)`.
- Lines 62-63: Absorption part of F3 signal; implemented by `f3_pos=f3_pos_pos+conj(f3_neg_neg)`.
- Lines 66-67: F2 Fourier transform; implemented by `f3f2_pos=fftshift(fft(f3_pos,parameters.zerofill(2),2),2)`.
- Lines 70-71: Absorption part of F2 signal; implemented by `f3f2=f3f2_pos+conj(f3f2_neg)`.
- Lines 73-74: F1 Fourier transform; implemented by `spectrum=fftshift(fft(f3f2,parameters.zerofill(1),1),1)`.
- Lines 76-78: Plotting; implemented by `kfigure(); plot_3d(spin_system,imag(spectrum),parameters, 10,[0.05 0.25 0.05 0.25],2,'positive')`.

### Key state/data transformations

- Lines 11: computes `sys.isotopes` using `sys.isotopes={'13C','13C','13C','1H','1H'}`.
- Lines 12: computes `sys.labels` using `sys.labels={'C','CA','CB','HA','HB'}`.
- Lines 15: computes `sys.magnet` using `sys.magnet=14.1`.
- Lines 16: computes `inter.zeeman.scalar` using `inter.zeeman.scalar={180 60 40 4 2}`.
- Lines 17: computes `inter.coupling.scalar{2,3}` using `inter.coupling.scalar{2,3}=35.0`.
- Lines 18: computes `inter.coupling.scalar{1,3}` using `inter.coupling.scalar{1,3}=0.3`.
- Lines 19: computes `inter.coupling.scalar{1,2}` using `inter.coupling.scalar{1,2}=55.0`.
- Lines 20: computes `inter.coupling.scalar{3,5}` using `inter.coupling.scalar{3,5}=130.0`.
- Lines 21: computes `inter.coupling.scalar{4,5}` using `inter.coupling.scalar{4,5}=7.0`.
- Lines 22: computes `inter.coupling.scalar{3,4}` using `inter.coupling.scalar{3,4}=6.0`.
- Lines 23: computes `inter.coupling.scalar{2,4}` using `inter.coupling.scalar{2,4}=140.0`.
- Lines 24: computes `inter.coupling.scalar{1,4}` using `inter.coupling.scalar{1,4}=4.0`.
- Lines 25: computes `inter.coupling.scalar{5,5}` using `inter.coupling.scalar{5,5}=0`.
- Lines 28: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 29: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 32: computes `parameters.J_ch` using `parameters.J_ch=140`.
- Lines 33: computes `parameters.J_cc` using `parameters.J_cc=35`.
- Lines 34: computes `parameters.delta` using `parameters.delta=1.1e-3`.

## Implementation structure

- 3D HCCH COSY experiment on a small protein fragment.
- Calculation time: seconds.
- Spin system
- Interactions
- Basis set
- Sequence parameters
- Spinach housekeeping
- Simulation
- Apodisation
- F3 Fourier transform
- Absorption part of F3 signal
- F2 Fourier transform

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `liquid()`, `apodisation()`, `fftshift()`, `conj()`, `kfigure()`, `plot_3d()`.
