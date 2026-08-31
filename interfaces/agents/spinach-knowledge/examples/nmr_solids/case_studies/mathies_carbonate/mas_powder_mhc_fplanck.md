# examples/nmr_solids/case_studies/mathies_carbonate/mas_powder_mhc_fplanck.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_solids/case_studies/mathies_carbonate/mas_powder_mhc_fplanck.m`
- Signature: `mas_powder_mhc_fplanck()`
- Total lines: 88

## Purpose

All protons in the unit cell of monohydrocalcite, magic angle spinning NMR simulation. Further details in: Calculation time: hours, minutes on a GPU.

## Physical / mathematical content

- Solid-state NMR examples. The key physics is anisotropic spin interactions under static or magic-angle-spinning conditions: chemical-shift anisotropy, dipolar coupling, quadrupolar coupling, cross-polarisation, and orientation averaging using Floquet, Fokker-Planck, or direct powder quadrature formalisms.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 12-13: 400 MHz NMR; implemented by `sys.magnet=9.4`.
- Lines 15-16: Read CASTEP file; implemented by `props=c2spinach('mhc.magres')`.
- Lines 18-19: Drop C, O, and Ca atoms; implemented by `drop_mask=ismember(props.symbols,{'C','O','Ca'})`.
- Lines 24-25: Get isotopes; implemented by `sys.isotopes={}`.
- Lines 34-36: Convert shielding tensors into shift using the parametrisation of Huang et al. ACIE 2021; implemented by `inter.zeeman.matrix=cell(1,numel(props.cst))`.
- Lines 43-44: Get coordinates; implemented by `inter.coordinates=mat2cell(props.std_geom,ones(18,1))`.
- Lines 46-47: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 52-53: Interaction cut-off, Hz; implemented by `sys.tols.inter_cutoff=500`.
- Lines 55-59: This needs a GPU sys.enable={'gpu'};; implemented by `spin_system=create(sys,inter)`.
- Lines 58-59: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 62-63: Pulse-acquire setup; implemented by `parameters.rate=10000`.
- Lines 75-76: Simulation; implemented by `fid=singlerot(spin_system,@acquire,parameters,'nmr')`.
- Lines 78-79: Apodisation; implemented by `fid=apodisation(spin_system,fid,{{'exp',6}})`.
- Lines 81-82: Fourier transform; implemented by `spectrum=fftshift(fft(fid,parameters.zerofill))`.
- Lines 84-85: Plotting; implemented by `kfigure(); plot_1d(spin_system,real(spectrum),parameters)`.

### Control flow inferred from the code

- Line 26: `for` loop over `n=1:numel(props.symbols)`.
- Line 27: conditional branch on `strcmp(props.symbols{n},'H')`.
- Line 37: `for` loop over `n=1:numel(props.cst)`.
- Line 38: conditional branch on `strcmp(sys.isotopes{n},'1H')`.

### Key state/data transformations

- Lines 13: computes `sys.magnet` using `sys.magnet=9.4`.
- Lines 16: computes `props` using `props=c2spinach('mhc.magres')`.
- Lines 19: computes `drop_mask` using `drop_mask=ismember(props.symbols,{'C','O','Ca'})`.
- Lines 20: computes `props.symbols(drop_mask)` using `props.symbols(drop_mask)=[]`.
- Lines 21: computes `props.std_geom(drop_mask,:)` using `props.std_geom(drop_mask,:)=[]`.
- Lines 22: computes `props.cst(drop_mask)` using `props.cst(drop_mask)=[]`.
- Lines 25: computes `sys.isotopes` using `sys.isotopes={}`.
- Lines 28: computes `sys.isotopes{n}` using `sys.isotopes{n}='1H'`.
- Lines 36: computes `inter.zeeman.matrix` using `inter.zeeman.matrix=cell(1,numel(props.cst))`.
- Lines 39: computes `inter.zeeman.matrix{n}` using `inter.zeeman.matrix{n}=29.25*eye(3)-props.cst{n}`.
- Lines 44: computes `inter.coordinates` using `inter.coordinates=mat2cell(props.std_geom,ones(18,1))`.
- Lines 47: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 48: computes `bas.approximation` using `bas.approximation='IK-0'`.
- Lines 49: computes `bas.level` using `bas.level=3`.
- Lines 50: computes `bas.projections` using `bas.projections=+1`.
- Lines 53: computes `sys.tols.inter_cutoff` using `sys.tols.inter_cutoff=500`.
- Lines 59: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 63: computes `parameters.rate` using `parameters.rate=10000`.

## Implementation structure

- All protons in the unit cell of monohydrocalcite, magic
- angle spinning NMR simulation. Further details in:
- Calculation time: hours, minutes on a GPU.
- 400 MHz NMR
- Read CASTEP file
- Drop C, O, and Ca atoms
- Get isotopes
- Convert shielding tensors into shift using the
- parametrisation of Huang et al. ACIE 2021
- Get coordinates
- Basis set
- Interaction cut-off, Hz

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `c2spinach()`, `ismember()`, `strcmp()`, `mat2cell()`, `create()`, `basis()`, `state()`, `singlerot()`, `apodisation()`, `fftshift()`, `kfigure()`, `plot_1d()`.
