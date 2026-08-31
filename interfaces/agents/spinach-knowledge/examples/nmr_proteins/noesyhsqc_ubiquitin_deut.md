# examples/nmr_proteins/noesyhsqc_ubiquitin_deut.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_proteins/noesyhsqc_ubiquitin_deut.m`
- Signature: `noesyhsqc_ubiquitin_deut()`
- Total lines: 101

## Purpose

1H-1H-15N NOESY-HSQC spectrum of 15N-labelled ubiquitin at 900 MHz with 90 ms mixing time. It is assumed that the protein is not 13C-labelled. Specific positions are deuterated, and deu- terium nuclei are simulated explicitly as spin-1 particles. Calculation time: a week on 32 cores, needs 512 GB of RAM.

## Physical / mathematical content

- Protein NMR examples. These files specialise liquid-state pulse sequences to labelled biomolecules, exploiting one-bond and two-bond heteronuclear couplings, coherence pathway filtering, selective decoupling, and high-dimensional indirect detection.
- The relaxation model is Redfield-type perturbation theory: fluctuating interactions enter through correlation functions or spectral densities and generate a linear relaxation superoperator.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 14-15: Protein data import; implemented by `options.pdb_mol=1`.
- Lines 23-24: Magnet field; implemented by `sys.magnet=21.1356`.
- Lines 26-27: Tolerances; implemented by `sys.tols.inter_cutoff=2.0`.
- Lines 30-31: Relaxation theory; implemented by `inter.relaxation={'redfield','t1_t2'}`.
- Lines 38-39: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 44-45: Algorithmic options; implemented by `sys.enable={'prop_cache','greedy'}`.
- Lines 48-49: Create the spin system structure; implemented by `spin_system=create(sys,inter)`.
- Lines 51-52: Kill carbons; implemented by `spin_system=kill_spin(spin_system,strcmp('13C',spin_system.comp.isotopes))`.
- Lines 54-55: Build the basis; implemented by `spin_system=basis(spin_system,bas)`.
- Lines 57-58: Sequence parameters; implemented by `parameters.tmix=0.090`.
- Lines 67-68: Simulation; implemented by `fid=liquid(spin_system,@noesyhsqc,parameters,'nmr')`.
- Lines 70-71: Apodisation; implemented by `fid.pos_pos=apodisation(spin_system,fid.pos_pos,{{'sqcos'},{'sqcos'},{'sqcos'}})`.
- Lines 76-77: F3 Fourier transform; implemented by `f3_pos_pos=fftshift(fft(fid.pos_pos,parameters.zerofill(3),3),3)`.
- Lines 82-83: Absorption part of F3 signal; implemented by `f3_pos=f3_pos_pos+conj(f3_neg_neg)`.
- Lines 86-87: F2 Fourier transform; implemented by `f3f2_pos=fftshift(fft(f3_pos,parameters.zerofill(2),2),2)`.
- Lines 90-91: Absorption part of F2 signal; implemented by `f3f2=f3f2_pos+conj(f3f2_neg)`.
- Lines 93-94: F1 Fourier transform; implemented by `spectrum=fftshift(fft(f3f2,parameters.zerofill(1),1),1)`.
- Lines 96-98: Plotting; implemented by `kfigure(); plot_3d(spin_system,-real(spectrum),parameters, 10,[0.01 0.1 0.01 0.1],2,'positive')`.

### Key state/data transformations

- Lines 15: computes `options.pdb_mol` using `options.pdb_mol=1`.
- Lines 16: computes `options.select` using `options.select='all'`.
- Lines 17: computes `options.noshift` using `options.noshift='delete'`.
- Lines 18-20: computes `options.deuterate` using `options.deuterate={'HA','HB','HB1','HB2','HB3','HG','HG1','HG2','HG3', 'HD','HD1','HD2','HD3','HE','HE1','HE2','HE3','HZ', 'HZ1','HZ2','HZ3','HH','HH1','HH2','HH3'}`.
- Lines 21: computes `[sys,inter]` using `[sys,inter]=protein('1D3Z.pdb','1D3Z.bmrb',options)`.
- Lines 24: computes `sys.magnet` using `sys.magnet=21.1356`.
- Lines 27: computes `sys.tols.inter_cutoff` using `sys.tols.inter_cutoff=2.0`.
- Lines 28: computes `sys.tols.prox_cutoff` using `sys.tols.prox_cutoff=4.0`.
- Lines 31: computes `inter.relaxation` using `inter.relaxation={'redfield','t1_t2'}`.
- Lines 32: computes `inter.r1_rates` using `inter.r1_rates=num2cell(100*strcmp('2H',sys.isotopes))`.
- Lines 33: computes `inter.r2_rates` using `inter.r2_rates=num2cell(100*strcmp('2H',sys.isotopes))`.
- Lines 34: computes `inter.rlx_keep` using `inter.rlx_keep='kite'`.
- Lines 35: computes `inter.equilibrium` using `inter.equilibrium='zero'`.
- Lines 36: computes `inter.tau_c` using `inter.tau_c={1e-8}`.
- Lines 39: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 40: computes `bas.approximation` using `bas.approximation='IK-1'`.
- Lines 41: computes `bas.connectivity` using `bas.connectivity='scalar_couplings'`.
- Lines 42: computes `bas.level` using `bas.level=4; bas.space_level=3`.

## Implementation structure

- 1H-1H-15N NOESY-HSQC spectrum of 15N-labelled ubiquitin at 900
- MHz with 90 ms mixing time. It is assumed that the protein is
- not 13C-labelled. Specific positions are deuterated, and deu-
- terium nuclei are simulated explicitly as spin-1 particles.
- Calculation time: a week on 32 cores, needs 512 GB of RAM.
- Protein data import
- Magnet field
- Tolerances
- Relaxation theory
- Basis set
- Algorithmic options
- Create the spin system structure

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `protein()`, `num2cell()`, `strcmp()`, `create()`, `kill_spin()`, `basis()`, `liquid()`, `apodisation()`, `fftshift()`, `conj()`, `kfigure()`, `plot_3d()`.
