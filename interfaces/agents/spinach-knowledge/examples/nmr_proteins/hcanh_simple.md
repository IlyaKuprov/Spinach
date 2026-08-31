# examples/nmr_proteins/hcanh_simple.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_proteins/hcanh_simple.m`
- Signature: `hcanh_simple()`
- Total lines: 81

## Purpose

A minimal example of H(CA)NH pulse sequence simulation. Calculation time: seconds.

## Physical / mathematical content

- Protein NMR examples. These files specialise liquid-state pulse sequences to labelled biomolecules, exploiting one-bond and two-bond heteronuclear couplings, coherence pathway filtering, selective decoupling, and high-dimensional indirect detection.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 10-11: Magnet field; implemented by `sys.magnet=14.1`.
- Lines 13-14: Spin system; implemented by `sys.isotopes={'15N','13C','1H','13C','1H'}`.
- Lines 17-18: Interactions; implemented by `inter.zeeman.scalar={110 60 8 180 4}`.
- Lines 31-32: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 35-36: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 39-40: Sequence parameters; implemented by `parameters.spins={'1H','15N','1H'}`.
- Lines 47-48: Simulation; implemented by `fid=liquid(spin_system,@hcanh,parameters,'nmr')`.
- Lines 50-51: Apodisation; implemented by `fid.pos_pos=apodisation(spin_system,fid.pos_pos,{{'sqcos'},{'sqcos'},{'sqcos'}})`.
- Lines 56-57: F3 Fourier transform; implemented by `f3_pos_pos=fftshift(fft(fid.pos_pos,parameters.zerofill(3),3),3)`.
- Lines 62-63: Absorption part of F3 signal; implemented by `f3_pos=f3_pos_pos+conj(f3_neg_neg)`.
- Lines 66-67: F2 Fourier transform; implemented by `f3f2_pos=fftshift(fft(f3_pos,parameters.zerofill(2),2),2)`.
- Lines 70-71: Absorption part of F2 signal; implemented by `f3f2=f3f2_pos+conj(f3f2_neg)`.
- Lines 73-74: F1 Fourier transform; implemented by `spectrum=fftshift(fft(f3f2,parameters.zerofill(1),1),1)`.
- Lines 76-78: Plotting; implemented by `kfigure(); plot_3d(spin_system,real(spectrum),parameters, 10,[0.05 0.5 0.05 0.5],2,'positive')`.

### Key state/data transformations

- Lines 11: computes `sys.magnet` using `sys.magnet=14.1`.
- Lines 14: computes `sys.isotopes` using `sys.isotopes={'15N','13C','1H','13C','1H'}`.
- Lines 15: computes `sys.labels` using `sys.labels={'N','CA','H','C','HA'}`.
- Lines 18: computes `inter.zeeman.scalar` using `inter.zeeman.scalar={110 60 8 180 4}`.
- Lines 19: computes `inter.coupling.scalar` using `inter.coupling.scalar=cell(5)`.
- Lines 20: computes `inter.coupling.scalar{1,3}` using `inter.coupling.scalar{1,3}=92`.
- Lines 21: computes `inter.coupling.scalar{1,2}` using `inter.coupling.scalar{1,2}=11`.
- Lines 22: computes `inter.coupling.scalar{1,4}` using `inter.coupling.scalar{1,4}=15`.
- Lines 23: computes `inter.coupling.scalar{1,5}` using `inter.coupling.scalar{1,5}=1`.
- Lines 24: computes `inter.coupling.scalar{2,4}` using `inter.coupling.scalar{2,4}=55`.
- Lines 25: computes `inter.coupling.scalar{2,3}` using `inter.coupling.scalar{2,3}=2`.
- Lines 26: computes `inter.coupling.scalar{2,5}` using `inter.coupling.scalar{2,5}=140`.
- Lines 27: computes `inter.coupling.scalar{3,4}` using `inter.coupling.scalar{3,4}=4`.
- Lines 28: computes `inter.coupling.scalar{3,5}` using `inter.coupling.scalar{3,5}=8`.
- Lines 29: computes `inter.coupling.scalar{4,5}` using `inter.coupling.scalar{4,5}=4`.
- Lines 32: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 33: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 36: computes `spin_system` using `spin_system=create(sys,inter)`.

## Implementation structure

- A minimal example of H(CA)NH pulse sequence simulation.
- Calculation time: seconds.
- Magnet field
- Spin system
- Interactions
- Basis set
- Spinach housekeeping
- Sequence parameters
- Simulation
- Apodisation
- F3 Fourier transform
- Absorption part of F3 signal

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `liquid()`, `apodisation()`, `fftshift()`, `conj()`, `kfigure()`, `plot_3d()`.
