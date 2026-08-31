# examples/nmr_solids/case_studies/mathies_14n_13c/mas_powder_gly_14n.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_solids/case_studies/mathies_14n_13c/mas_powder_gly_14n.m`
- Signature: `mas_powder_gly_14n()`
- Total lines: 91

## Purpose

14N MAS spectrum of glycine powder (assuming decoupling of 1H and 13C), computed using the Fokker-Planck MAS formalism and a spherical grid. Numerical rotating frame transformation is used because 14N quadrupolar interaction is large. Calculation time: seconds.

## Physical / mathematical content

- Solid-state NMR examples. The key physics is anisotropic spin interactions under static or magic-angle-spinning conditions: chemical-shift anisotropy, dipolar coupling, quadrupolar coupling, cross-polarisation, and orientation averaging using Floquet, Fokker-Planck, or direct powder quadrature formalisms.
- The file uses a Fokker-Planck-style enlarged state space in which spatial or orientational coordinates are promoted to extra dimensions and coupled to spin dynamics through differential operators.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.
- Quadrupolar physics is relevant: nuclei with spin > 1/2 interact with the electric field gradient tensor, introducing second-rank anisotropy, asymmetry, and overtone or MQ phenomena.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 14-15: Read CASTEP file; implemented by `props=c2spinach('glycine.magres')`.
- Lines 17-18: Drop H, O, and C atoms; implemented by `drop_mask=ismember(props.symbols,{'H','O','C'})`.
- Lines 24-25: keep only 1 14N; implemented by `sys.isotopes{1}='14N'`.
- Lines 27-28: Convert shielding tensor into shift; implemented by `inter.zeeman.matrix{1}=-props.cst{1}`.
- Lines 30-31: Set isotropic chemical shift to experimental value; implemented by `inter.zeeman.matrix=shift_iso(inter.zeeman.matrix,1,110.0)`.
- Lines 33-34: Quadrupolar interaction from CASTEP; implemented by `nqi=castep2nqi(props.efg{1},20.44e-3,1)`.
- Lines 37-38: Magnet field; implemented by `sys.magnet=9.4`.
- Lines 40-41: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 44-45: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 48-49: Experiment setup; implemented by `parameters.sweep=3e6`.
- Lines 60-61: Simulation; implemented by `fid=powder(spin_system,@acquire,parameters,'nmr')`.
- Lines 63-64: Apodisation; implemented by `fid=apodisation(spin_system,fid,{{'exp',6}})`.
- Lines 66-67: Fourier transform; implemented by `spectrum=fftshift(fft(fid,parameters.zerofill))`.
- Lines 69-70: Plotting; implemented by `kfigure(); hold on`.
- Lines 73-74: Quadrupolar interaction measured by O'Dell PCCP 2009; implemented by `dell_nqi=2*pi*eeqq2nqi(1.18e6,0.53,1,[0 0 0])`.
- Lines 86-87: Plotting and legend; implemented by `plot_1d(spin_system,real(spectrum),parameters)`.

### Key state/data transformations

- Lines 15: computes `props` using `props=c2spinach('glycine.magres')`.
- Lines 18: computes `drop_mask` using `drop_mask=ismember(props.symbols,{'H','O','C'})`.
- Lines 19: computes `props.symbols(drop_mask)` using `props.symbols(drop_mask)=[]`.
- Lines 20: computes `props.std_geom(drop_mask,:)` using `props.std_geom(drop_mask,:)=[]`.
- Lines 21: computes `props.cst(drop_mask)` using `props.cst(drop_mask)=[]`.
- Lines 22: computes `props.efg(drop_mask)` using `props.efg(drop_mask)=[]`.
- Lines 25: computes `sys.isotopes{1}` using `sys.isotopes{1}='14N'`.
- Lines 28: computes `inter.zeeman.matrix{1}` using `inter.zeeman.matrix{1}=-props.cst{1}`.
- Lines 31: computes `inter.zeeman.matrix` using `inter.zeeman.matrix=shift_iso(inter.zeeman.matrix,1,110.0)`.
- Lines 34: computes `nqi` using `nqi=castep2nqi(props.efg{1},20.44e-3,1)`.
- Lines 35: computes `inter.coupling.matrix{1,1}` using `inter.coupling.matrix{1,1}=remtrace(nqi)`.
- Lines 38: computes `sys.magnet` using `sys.magnet=9.4`.
- Lines 41: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 42: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 45: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 49: computes `parameters.sweep` using `parameters.sweep=3e6`.
- Lines 50: computes `parameters.npoints` using `parameters.npoints=256`.
- Lines 51: computes `parameters.zerofill` using `parameters.zerofill=1024`.

## Implementation structure

- 14N MAS spectrum of glycine powder (assuming decoupling of 1H
- and 13C), computed using the Fokker-Planck MAS formalism and
- a spherical grid. Numerical rotating frame transformation is
- used because 14N quadrupolar interaction is large.
- Calculation time: seconds.
- Read CASTEP file
- Drop H, O, and C atoms
- keep only 1 14N
- Convert shielding tensor into shift
- Set isotropic chemical shift to experimental value
- Quadrupolar interaction from CASTEP
- Magnet field

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `c2spinach()`, `ismember()`, `shift_iso()`, `castep2nqi()`, `remtrace()`, `create()`, `basis()`, `state()`, `powder()`, `apodisation()`, `fftshift()`, `kfigure()`, `plot_1d()`, `eeqq2nqi()`, `klegend()`.
