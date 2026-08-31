# examples/nmr_solids/case_studies/mathies_14n_13c/mas_powder_gly_13c.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_solids/case_studies/mathies_14n_13c/mas_powder_gly_13c.m`
- Signature: `mas_powder_gly_13c()`
- Total lines: 103

## Purpose

13C MAS spectrum of glycine powder (assuming decoupling of 1H), computed using the Fokker-Planck MAS formalism and a spherical grid. The field dependence of the line shape of 13CA due to the presence of the quadrupolar 14N nucleus is shown. The calcula- tion is performed in the rotating frame with respect to 13C and the laboratory frame with respect to 14N. Calculation time: seconds.

## Physical / mathematical content

- Solid-state NMR examples. The key physics is anisotropic spin interactions under static or magic-angle-spinning conditions: chemical-shift anisotropy, dipolar coupling, quadrupolar coupling, cross-polarisation, and orientation averaging using Floquet, Fokker-Planck, or direct powder quadrature formalisms.
- The file uses a Fokker-Planck-style enlarged state space in which spatial or orientational coordinates are promoted to extra dimensions and coupled to spin dynamics through differential operators.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.
- Quadrupolar physics is relevant: nuclei with spin > 1/2 interact with the electric field gradient tensor, introducing second-rank anisotropy, asymmetry, and overtone or MQ phenomena.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 16-17: Read CASTEP file; implemented by `props=c2spinach('glycine.magres')`.
- Lines 19-20: Drop H and O atoms; implemented by `drop_mask=ismember(props.symbols,{'H','O'})`.
- Lines 26-27: Keep 13CA and 14N; implemented by `sys.isotopes{1}='13C'`.
- Lines 30-31: Convert shielding tensors into shift; implemented by `inter.zeeman.matrix{1}=-props.cst{2}`.
- Lines 34-35: Set isotropic chemical shifts to experimental values; implemented by `inter.zeeman.matrix=shift_iso(inter.zeeman.matrix,[1 2],[43.6 110.0])`.
- Lines 37-38: Cartesian coordinates; implemented by `inter.coordinates={props.std_geom(2,:)`.
- Lines 41-42: Quadrupolar interaction; implemented by `inter.coupling.matrix=cell(2,2)`.
- Lines 46-47: Experiment setup; implemented by `parameters.rate=10000`.
- Lines 56-57: Numerical rotating frame transforms; implemented by `parameters.rframes={{'13C',1},{'14N',2}}`.
- Lines 59-60: Vary the field; implemented by `fields=[4.7 9.4 14.1]`.
- Lines 62-63: Get a figure going; implemented by `kfigure(); hold on`.
- Lines 67-68: Magnet field; implemented by `sys.magnet=fields(m)`.
- Lines 70-71: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 74-75: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 78-79: Experiment setup; implemented by `parameters.sweep=200*sys.magnet`.
- Lines 84-85: Lab frame Hamiltonian, then numerical rotating frames; implemented by `fid=singlerot(spin_system,@acquire,parameters,'labframe')`.
- Lines 87-88: Apodisation; implemented by `fid=apodisation(spin_system,fid,{{'exp',6}})`.
- Lines 90-91: Fourier transform; implemented by `spectrum=fftshift(fft(fid,parameters.zerofill))`.

### Control flow inferred from the code

- Line 65: `for` loop over `m=1:numel(fields)`.

### Key state/data transformations

- Lines 17: computes `props` using `props=c2spinach('glycine.magres')`.
- Lines 20: computes `drop_mask` using `drop_mask=ismember(props.symbols,{'H','O'})`.
- Lines 21: computes `props.symbols(drop_mask)` using `props.symbols(drop_mask)=[]`.
- Lines 22: computes `props.std_geom(drop_mask,:)` using `props.std_geom(drop_mask,:)=[]`.
- Lines 23: computes `props.cst(drop_mask)` using `props.cst(drop_mask)=[]`.
- Lines 24: computes `props.efg(drop_mask)` using `props.efg(drop_mask)=[]`.
- Lines 27: computes `sys.isotopes{1}` using `sys.isotopes{1}='13C'`.
- Lines 28: computes `sys.isotopes{2}` using `sys.isotopes{2}='14N'`.
- Lines 31: computes `inter.zeeman.matrix{1}` using `inter.zeeman.matrix{1}=-props.cst{2}`.
- Lines 32: computes `inter.zeeman.matrix{2}` using `inter.zeeman.matrix{2}=-props.cst{5}`.
- Lines 35: computes `inter.zeeman.matrix` using `inter.zeeman.matrix=shift_iso(inter.zeeman.matrix,[1 2],[43.6 110.0])`.
- Lines 38: computes `inter.coordinates` using `inter.coordinates={props.std_geom(2,:)`.
- Lines 42: computes `inter.coupling.matrix` using `inter.coupling.matrix=cell(2,2)`.
- Lines 43: computes `nqi` using `nqi=castep2nqi(props.efg{5},20.44e-3,1)`.
- Lines 44: computes `inter.coupling.matrix{2,2}` using `inter.coupling.matrix{2,2}=remtrace(nqi)`.
- Lines 47: computes `parameters.rate` using `parameters.rate=10000`.
- Lines 48: computes `parameters.axis` using `parameters.axis=[1 1 1]`.
- Lines 49: computes `parameters.max_rank` using `parameters.max_rank=5`.

## Implementation structure

- 13C MAS spectrum of glycine powder (assuming decoupling of 1H),
- computed using the Fokker-Planck MAS formalism and a spherical
- grid. The field dependence of the line shape of 13CA due to the
- presence of the quadrupolar 14N nucleus is shown. The calcula-
- tion is performed in the rotating frame with respect to 13C and
- the laboratory frame with respect to 14N.
- Calculation time: seconds.
- Read CASTEP file
- Drop H and O atoms
- Keep 13CA and 14N
- Convert shielding tensors into shift
- Set isotropic chemical shifts to experimental values

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `c2spinach()`, `ismember()`, `shift_iso()`, `castep2nqi()`, `remtrace()`, `kfigure()`, `fields()`, `create()`, `basis()`, `spin()`, `state()`, `singlerot()`, `apodisation()`, `fftshift()`, `plot_1d()`, `klegend()`.
