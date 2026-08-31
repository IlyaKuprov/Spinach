# examples/nmr_solids/case_studies/mathies_carbonate/mas_powder_mhc_fplanck_exchange.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_solids/case_studies/mathies_carbonate/mas_powder_mhc_fplanck_exchange.m`
- Signature: `mas_powder_mhc_fplanck_exchange()`
- Total lines: 87

## Purpose

Water protons in the unit cell of monohydrocalcite, inc- luding position exchange and MAS. Further details in: Calculation time: seconds.

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
- Lines 24-26: Two reaction endpoints with two protons each, swapped by the reaction; implemented by `sys.isotopes{1}='1H'`.
- Lines 31-33: Convert shielding tensors into shift using the parametrisation of Huang et al. ACIE 2021; implemented by `inter.zeeman.matrix{1}=29.25*eye(3)-props.cst{1}`.
- Lines 38-39: Get coordinates; implemented by `inter.coordinates{1}=props.std_geom(1,:)`.
- Lines 44-45: Chemical kinetics endpoints; implemented by `inter.chem.parts={[1 2],[3 4]}`.
- Lines 47-48: Reaction rate matrix, Hz; implemented by `inter.chem.rates=2e3*[-1 1`.
- Lines 51-52: Initial concentrations (arb. units); implemented by `inter.chem.concs=[1 1]`.
- Lines 54-55: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 58-59: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 62-63: Pulse-acquire setup; implemented by `parameters.rate=10000`.
- Lines 74-75: Simulation; implemented by `fid=singlerot(spin_system,@acquire,parameters,'nmr')`.
- Lines 77-78: Apodisation; implemented by `fid=apodisation(spin_system,fid,{{'exp',6}})`.
- Lines 80-81: Fourier transform; implemented by `spectrum=fftshift(fft(fid,parameters.zerofill))`.
- Lines 83-84: Plotting; implemented by `kfigure(); plot_1d(spin_system,real(spectrum),parameters)`.

### Key state/data transformations

- Lines 13: computes `sys.magnet` using `sys.magnet=9.4`.
- Lines 16: computes `props` using `props=c2spinach('mhc.magres')`.
- Lines 19: computes `drop_mask` using `drop_mask=ismember(props.symbols,{'C','O','Ca'})`.
- Lines 20: computes `props.symbols(drop_mask)` using `props.symbols(drop_mask)=[]`.
- Lines 21: computes `props.std_geom(drop_mask,:)` using `props.std_geom(drop_mask,:)=[]`.
- Lines 22: computes `props.cst(drop_mask)` using `props.cst(drop_mask)=[]`.
- Lines 26: computes `sys.isotopes{1}` using `sys.isotopes{1}='1H'`.
- Lines 27: computes `sys.isotopes{2}` using `sys.isotopes{2}='1H'`.
- Lines 28: computes `sys.isotopes{3}` using `sys.isotopes{3}='1H'`.
- Lines 29: computes `sys.isotopes{4}` using `sys.isotopes{4}='1H'`.
- Lines 33: computes `inter.zeeman.matrix{1}` using `inter.zeeman.matrix{1}=29.25*eye(3)-props.cst{1}`.
- Lines 34: computes `inter.zeeman.matrix{2}` using `inter.zeeman.matrix{2}=29.25*eye(3)-props.cst{4}`.
- Lines 35: computes `inter.zeeman.matrix{3}` using `inter.zeeman.matrix{3}=29.25*eye(3)-props.cst{4}`.
- Lines 36: computes `inter.zeeman.matrix{4}` using `inter.zeeman.matrix{4}=29.25*eye(3)-props.cst{1}`.
- Lines 39: computes `inter.coordinates{1}` using `inter.coordinates{1}=props.std_geom(1,:)`.
- Lines 40: computes `inter.coordinates{2}` using `inter.coordinates{2}=props.std_geom(4,:)`.
- Lines 41: computes `inter.coordinates{3}` using `inter.coordinates{3}=props.std_geom(4,:)`.
- Lines 42: computes `inter.coordinates{4}` using `inter.coordinates{4}=props.std_geom(1,:)`.

## Implementation structure

- Water protons in the unit cell of monohydrocalcite, inc-
- luding position exchange and MAS. Further details in:
- Calculation time: seconds.
- 400 MHz NMR
- Read CASTEP file
- Drop C, O, and Ca atoms
- Two reaction endpoints with two protons
- each, swapped by the reaction
- Convert shielding tensors into shift using the
- parametrisation of Huang et al. ACIE 2021
- Get coordinates
- Chemical kinetics endpoints

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `c2spinach()`, `ismember()`, `create()`, `basis()`, `state()`, `singlerot()`, `apodisation()`, `fftshift()`, `kfigure()`, `plot_1d()`.
