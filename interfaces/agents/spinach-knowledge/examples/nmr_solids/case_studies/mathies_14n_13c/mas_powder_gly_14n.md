# examples/nmr_solids/case_studies/mathies_14n_13c/mas_powder_gly_14n.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_solids/case_studies/mathies_14n_13c/mas_powder_gly_14n.m`
- Signature: `mas_powder_gly_14n()`
- Total lines: 93

## Purpose

Static 14N powder spectrum of glycine (assuming decoupling of 1H and 13C), computed on a spherical grid from CASTEP tensors and compared to a simulation using the quadrupolar parameters measured by O'Dell and Schurko (PCCP 2009, DOI 10.1039/b906114b). Numerical rotating frame transformation is used because 14N quadrupolar interaction is large. Calculation time: seconds.

## Physical / mathematical content

- Solid-state NMR examples. The key physics is anisotropic spin interactions under static conditions: chemical-shift anisotropy, quadrupolar coupling, and orientation averaging using direct powder quadrature on a spherical grid.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.
- Quadrupolar physics is relevant: nuclei with spin > 1/2 interact with the electric field gradient tensor, introducing second-rank anisotropy, asymmetry, and overtone or MQ phenomena.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 16-17: Read CASTEP file; implemented by `props=c2spinach('glycine.magres')`.
- Lines 19-20: Drop H, O, and C atoms; implemented by `drop_mask=ismember(props.symbols,{'H','O','C'})`.
- Lines 26-27: keep only 1 14N; implemented by `sys.isotopes{1}='14N'`.
- Lines 29-30: Convert shielding tensor into shift; implemented by `inter.zeeman.matrix{1}=-props.cst{1}`.
- Lines 32-33: Set isotropic chemical shift to experimental value; implemented by `inter.zeeman.matrix=shift_iso(inter.zeeman.matrix,1,110.0)`.
- Lines 35-36: Quadrupolar interaction from CASTEP; implemented by `nqi=castep2nqi(props.efg{1},20.44e-3,1)`.
- Lines 39-40: Magnet field; implemented by `sys.magnet=9.4`.
- Lines 42-43: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 46-47: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 50-51: Experiment setup; implemented by `parameters.sweep=3e6`.
- Lines 62-63: Simulation; implemented by `fid=powder(spin_system,@acquire,parameters,'nmr')`.
- Lines 65-66: Apodisation; implemented by `fid=apodisation(spin_system,fid,{{'exp',6}})`.
- Lines 68-69: Fourier transform; implemented by `spectrum=fftshift(fft(fid,parameters.zerofill))`.
- Lines 71-72: Plotting; implemented by `kfigure(); hold on`.
- Lines 75-76: Quadrupolar interaction measured by O'Dell PCCP 2009; implemented by `dell_nqi=2*pi*eeqq2nqi(1.18e6,0.53,1,[0 0 0])`.
- Lines 88-89: Plotting and legend; implemented by `plot_1d(spin_system,real(spectrum),parameters)`.

### Key state/data transformations

- Lines 17: computes `props` using `props=c2spinach('glycine.magres')`.
- Lines 20: computes `drop_mask` using `drop_mask=ismember(props.symbols,{'H','O','C'})`.
- Lines 21: computes `props.symbols(drop_mask)` using `props.symbols(drop_mask)=[]`.
- Lines 22: computes `props.std_geom(drop_mask,:)` using `props.std_geom(drop_mask,:)=[]`.
- Lines 23: computes `props.cst(drop_mask)` using `props.cst(drop_mask)=[]`.
- Lines 24: computes `props.efg(drop_mask)` using `props.efg(drop_mask)=[]`.
- Lines 27: computes `sys.isotopes{1}` using `sys.isotopes{1}='14N'`.
- Lines 30: computes `inter.zeeman.matrix{1}` using `inter.zeeman.matrix{1}=-props.cst{1}`.
- Lines 33: computes `inter.zeeman.matrix` using `inter.zeeman.matrix=shift_iso(inter.zeeman.matrix,1,110.0)`.
- Lines 36: computes `nqi` using `nqi=castep2nqi(props.efg{1},20.44e-3,1)`.
- Lines 37: computes `inter.coupling.matrix{1,1}` using `inter.coupling.matrix{1,1}=remtrace(nqi)`.
- Lines 40: computes `sys.magnet` using `sys.magnet=9.4`.
- Lines 43: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 44: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 47: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 51: computes `parameters.sweep` using `parameters.sweep=3e6`.
- Lines 52: computes `parameters.npoints` using `parameters.npoints=256`.
- Lines 53: computes `parameters.zerofill` using `parameters.zerofill=1024`.

## Implementation structure

- Static 14N powder spectrum of glycine (assuming decoupling of
- 1H and 13C), computed on a spherical grid from CASTEP tensors
- and compared to a simulation using the quadrupolar parameters
- measured by O'Dell and Schurko (PCCP 2009, DOI 10.1039/b906114b).
- Numerical rotating frame transformation is used because 14N
- quadrupolar interaction is large.
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
