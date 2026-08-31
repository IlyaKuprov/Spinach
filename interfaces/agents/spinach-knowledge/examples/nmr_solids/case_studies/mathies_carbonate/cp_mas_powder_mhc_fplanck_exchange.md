# examples/nmr_solids/case_studies/mathies_carbonate/cp_mas_powder_mhc_fplanck_exchange.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_solids/case_studies/mathies_carbonate/cp_mas_powder_mhc_fplanck_exchange.m`
- Signature: `cp_mas_powder_mhc_fplanck_exchange()`
- Total lines: 115

## Purpose

Cross-polarisation contact curve under magic angle spinning in the presence of chemical exchange for H1, H4 and C19 in the unit cell of monohydrocalcite. Further details in: Calculation time: hours, much faster on a GPU.

## Physical / mathematical content

- Solid-state NMR examples. The key physics is anisotropic spin interactions under static or magic-angle-spinning conditions: chemical-shift anisotropy, dipolar coupling, quadrupolar coupling, cross-polarisation, and orientation averaging using Floquet, Fokker-Planck, or direct powder quadrature formalisms.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 13-14: 400 MHz NMR; implemented by `sys.magnet=9.4`.
- Lines 16-17: Read CASTEP file; implemented by `props=c2spinach('mhc.magres')`.
- Lines 19-20: Drop O and Ca atoms; implemented by `drop_mask=ismember(props.symbols,{'O','Ca'})`.
- Lines 25-27: Two chemical endpoints: H1, H4, and C19, with H1 and H4 under chemical exchange; implemented by `sys.isotopes{1}='1H'; sys.isotopes{4}='1H'`.
- Lines 31-33: Convert shielding tensors into shift using the parametrisation of Huang et al. ACIE 2021; implemented by `inter.zeeman.matrix{1}=29.25*eye(3)-props.cst{1}`.
- Lines 40-41: Cartesian coordinates; implemented by `inter.coordinates={props.std_geom(1,:)`.
- Lines 48-49: Chemical kinetics endpoints; implemented by `inter.chem.parts={[1 2 3],[4 5 6]}`.
- Lines 51-52: Initial concentrations (arb. units); implemented by `inter.chem.concs=[1 1]`.
- Lines 54-55: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 58-59: Disable start-up checks; implemented by `sys.disable={'hygiene'}`.
- Lines 61-65: Enable GPU sys.enable={'gpu'};; implemented by `exch_rates=[1e1 1e2 1e3 1e4 1e5 1e6]`.
- Lines 64-65: Exchange rate constant array; implemented by `exch_rates=[1e1 1e2 1e3 1e4 1e5 1e6]`.
- Lines 67-68: Experiment setup; implemented by `parameters.spins={'1H','13C'}`.
- Lines 81-83: Preallocate contact curve array; implemented by `contact_curves=zeros(numel(exch_rates), parameters.nsteps+1)`.
- Lines 85-86: Loop over exchange rates; implemented by `for n=1:numel(exch_rates)`.
- Lines 88-89: Set the exchange rates; implemented by `inter.chem.rates=exch_rates(n)*[-1 1`.
- Lines 92-93: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 96-97: Detection state; implemented by `parameters.coil=state(spin_system,'L+','13C')`.

### Control flow inferred from the code

- Line 86: `for` loop over `n=1:numel(exch_rates)`.

### Key state/data transformations

- Lines 14: computes `sys.magnet` using `sys.magnet=9.4`.
- Lines 17: computes `props` using `props=c2spinach('mhc.magres')`.
- Lines 20: computes `drop_mask` using `drop_mask=ismember(props.symbols,{'O','Ca'})`.
- Lines 21: computes `props.symbols(drop_mask)` using `props.symbols(drop_mask)=[]`.
- Lines 22: computes `props.std_geom(drop_mask,:)` using `props.std_geom(drop_mask,:)=[]`.
- Lines 23: computes `props.cst(drop_mask)` using `props.cst(drop_mask)=[]`.
- Lines 27: computes `sys.isotopes{1}` using `sys.isotopes{1}='1H'; sys.isotopes{4}='1H'`.
- Lines 28: computes `sys.isotopes{2}` using `sys.isotopes{2}='1H'; sys.isotopes{5}='1H'`.
- Lines 29: computes `sys.isotopes{3}` using `sys.isotopes{3}='13C'; sys.isotopes{6}='13C'`.
- Lines 33: computes `inter.zeeman.matrix{1}` using `inter.zeeman.matrix{1}=29.25*eye(3)-props.cst{1}`.
- Lines 34: computes `inter.zeeman.matrix{2}` using `inter.zeeman.matrix{2}=29.25*eye(3)-props.cst{4}`.
- Lines 35: computes `inter.zeeman.matrix{3}` using `inter.zeeman.matrix{3}=169.86*eye(3)-props.cst{19}`.
- Lines 36: computes `inter.zeeman.matrix{4}` using `inter.zeeman.matrix{4}=29.25*eye(3)-props.cst{4}`.
- Lines 37: computes `inter.zeeman.matrix{5}` using `inter.zeeman.matrix{5}=29.25*eye(3)-props.cst{1}`.
- Lines 38: computes `inter.zeeman.matrix{6}` using `inter.zeeman.matrix{6}=169.86*eye(3)-props.cst{19}`.
- Lines 41: computes `inter.coordinates` using `inter.coordinates={props.std_geom(1,:)`.
- Lines 49: computes `inter.chem.parts` using `inter.chem.parts={[1 2 3],[4 5 6]}`.
- Lines 52: computes `inter.chem.concs` using `inter.chem.concs=[1 1]`.

## Implementation structure

- Cross-polarisation contact curve under magic angle spinning
- in the presence of chemical exchange for H1, H4 and C19 in
- the unit cell of monohydrocalcite. Further details in:
- Calculation time: hours, much faster on a GPU.
- 400 MHz NMR
- Read CASTEP file
- Drop O and Ca atoms
- Two chemical endpoints: H1, H4, and C19,
- with H1 and H4 under chemical exchange
- Convert shielding tensors into shift using the
- parametrisation of Huang et al. ACIE 2021
- Cartesian coordinates

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `c2spinach()`, `ismember()`, `exch_rates()`, `create()`, `basis()`, `state()`, `contact_curves()`, `singlerot()`, `kfigure()`, `kxlabel()`, `kylabel()`.
