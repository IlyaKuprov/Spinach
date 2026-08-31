# examples/nmr_proteins/hnco_simple.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_proteins/hnco_simple.m`
- Signature: `hnco_simple()`
- Total lines: 77

## Purpose

A minimal example of HNCO pulse sequence simulation. Calculation time: seconds.

## Physical / mathematical content

- Protein NMR examples. These files specialise liquid-state pulse sequences to labelled biomolecules, exploiting one-bond and two-bond heteronuclear couplings, coherence pathway filtering, selective decoupling, and high-dimensional indirect detection.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 10-11: Magnet field; implemented by `sys.magnet=14.1`.
- Lines 13-14: Spin system; implemented by `sys.isotopes={'15N','13C','1H','13C'}`.
- Lines 17-18: Interactions; implemented by `inter.zeeman.scalar={110 55 8.0 180}`.
- Lines 25-26: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 29-30: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 33-34: Sequence parameters; implemented by `parameters.spins={'15N','13C','1H'}`.
- Lines 43-44: Simulation; implemented by `fid=liquid(spin_system,@hnco,parameters,'nmr')`.
- Lines 46-47: Apodisation; implemented by `fid.pos_pos=apodisation(spin_system,fid.pos_pos,{{'sqcos'},{'sqcos'},{'sqcos'}})`.
- Lines 52-53: F3 Fourier transform; implemented by `f3_pos_pos=fftshift(fft(fid.pos_pos,parameters.zerofill(3),3),3)`.
- Lines 58-59: Absorption part of F3 signal; implemented by `f3_pos=f3_pos_pos+conj(f3_neg_neg)`.
- Lines 62-63: F2 Fourier transform; implemented by `f3f2_pos=fftshift(fft(f3_pos,parameters.zerofill(2),2),2)`.
- Lines 66-67: Absorption part of F2 signal; implemented by `f3f2=f3f2_pos+conj(f3f2_neg)`.
- Lines 69-70: F1 Fourier transform; implemented by `spectrum=fftshift(fft(f3f2,parameters.zerofill(1),1),1)`.
- Lines 72-74: Plotting; implemented by `kfigure(); plot_3d(spin_system,real(spectrum),parameters, 10,[0.2 0.9 0.2 0.9],2,'positive')`.

### Key state/data transformations

- Lines 11: computes `sys.magnet` using `sys.magnet=14.1`.
- Lines 14: computes `sys.isotopes` using `sys.isotopes={'15N','13C','1H','13C'}`.
- Lines 15: computes `sys.labels` using `sys.labels={'N','CA','H','C'}`.
- Lines 18: computes `inter.zeeman.scalar` using `inter.zeeman.scalar={110 55 8.0 180}`.
- Lines 19: computes `inter.coupling.scalar` using `inter.coupling.scalar=cell(4)`.
- Lines 20: computes `inter.coupling.scalar{1,3}` using `inter.coupling.scalar{1,3}=92`.
- Lines 21: computes `inter.coupling.scalar{1,2}` using `inter.coupling.scalar{1,2}=11`.
- Lines 22: computes `inter.coupling.scalar{1,4}` using `inter.coupling.scalar{1,4}=15`.
- Lines 23: computes `inter.coupling.scalar{2,4}` using `inter.coupling.scalar{2,4}=55`.
- Lines 26: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 27: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 30: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 34: computes `parameters.spins` using `parameters.spins={'15N','13C','1H'}`.
- Lines 35: computes `parameters.offset` using `parameters.offset=[-7200 25000 4800]`.
- Lines 36: computes `parameters.sweep` using `parameters.sweep=[5000 10000 5000]`.
- Lines 37: computes `parameters.npoints` using `parameters.npoints=[63 64 65]`.
- Lines 38: computes `parameters.zerofill` using `parameters.zerofill=[255 256 257]`.
- Lines 39: computes `parameters.tau` using `parameters.tau=[2.25e-3 14e-3 4e-3]`.

## Implementation structure

- A minimal example of HNCO pulse sequence simulation.
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
