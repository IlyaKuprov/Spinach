# interfaces/g2spinach.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/interfaces/g2spinach.m`
- Signature: `[sys,inter]=g2spinach(props,particles,references,options)`
- Total lines: 290

## Purpose

Makes Spinach data structures from parsed outputs of electronic structure theory packages, such as Gaussian and ORCA. Syntax: [sys,inter]=g2spinach(props,particles,references,options)

## Physical / mathematical content

- This file belongs to the `interfaces` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.
- Quadrupolar physics is relevant: nuclei with spin > 1/2 interact with the electric field gradient tensor, introducing second-rank anisotropy, asymmetry, and overtone or MQ phenomena.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `numel()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 94-95: Set default options; implemented by `if ~exist('options','var'), options=[]; end`.
- Lines 97-98: Check consistency; implemented by `grumble(props,particles,references,options)`.
- Lines 100-101: Fundamental constants; implemented by `nuclear_magneton=7.6225932291E6`.
- Lines 103-104: Index the particles to include; implemented by `sys.isotopes={}; index=[]; ref_index=[]`.
- Lines 115-116: Decide whether to include coordinates; implemented by `if nargin<4`.
- Lines 128-129: Process coordinates; implemented by `if include_xyz`.
- Lines 134-135: Decide which magnetic parameters to return; implemented by `switch ismember('E',[particles{:}])`.
- Lines 137-138: EPR parameterization; implemented by `case 1`.
- Lines 140-141: Add the electron as the last spin in the isotope list; implemented by `sys.isotopes=[sys.isotopes, 'E']; nspins=nspins+1`.
- Lines 143-144: Electron coordinates should be treated as unknown; implemented by `if include_xyz, inter.coordinates{end+1}=[]; end`.
- Lines 146-147: All Zeeman tensors are zero except for electron g-tensor; implemented by `inter.zeeman.matrix=cell(1,nspins)`.
- Lines 152-153: All couplings are zero except for the hyperfine couplings to the electron; implemented by `inter.coupling.matrix=mat2cell(zeros(3*nspins,3*nspins),3*ones(nspins,1),3*ones(nspins,1))`.
- Lines 159-160: Remove small hyperfine couplings; implemented by `if (nargin==4)&&isfield(options,'min_hfc')`.
- Lines 170-171: Remove the nuclei that are not coupled to the electron; implemented by `if (nargin==4)&&isfield(options,'purge')&&strcmp(options.purge,'on')`.
- Lines 180-181: NMR parameterization; implemented by `case 0`.
- Lines 183-184: Reference to bare nuclei if the user did not supply any reference values; implemented by `if (nargin<3)||(~any(references)), references=zeros(size(particles)); end`.
- Lines 186-187: Absorb Zeeman tensors; implemented by `if isfield(props,'cst')`.
- Lines 194-195: Absorb quadrupolar couplings; implemented by `if isfield(props,'nqi')`.

### Control flow inferred from the code

- Line 95: conditional branch on `~exist('options','var'), options=[]; end`.
- Line 105: `for` loop over `n=1:length(props.symbols)`.
- Line 106: `for` loop over `k=1:length(particles)`.
- Line 107: conditional branch on `strcmp(props.symbols{n},particles{k}{1})`.
- Line 116: conditional branch on `nargin<4`.
- Line 119: conditional branch on `~isfield(options,'no_xyz')`.
- Line 129: conditional branch on `include_xyz`.
- Line 135: dispatches on `ismember('E',[particles{:}])`; cases `1`, `0`.
- Line 144: conditional branch on `include_xyz, inter.coordinates{end+1}=[]; end`.
- Line 154: `for` loop over `n=1:(nspins-1)`.
- Line 160: conditional branch on `(nargin==4)&&isfield(options,'min_hfc')`.
- Line 161: `for` loop over `n=1:nspins`.
- Line 162: `for` loop over `k=1:nspins`.
- Line 163: conditional branch on `norm(inter.coupling.matrix{n,k},'fro')<options.min_hfc`.

### Key state/data transformations

- Lines 101: computes `nuclear_magneton` using `nuclear_magneton=7.6225932291E6`.
- Lines 104: computes `sys.isotopes` using `sys.isotopes={}; index=[]; ref_index=[]`.
- Lines 109: computes `index` using `index=[index n]; ref_index=[ref_index k]`.
- Lines 113: computes `nspins` using `nspins=length(index)`.
- Lines 117: computes `include_xyz` using `include_xyz=true()`.
- Lines 130: computes `inter.coordinates` using `inter.coordinates=props.std_geom(index,:)`.
- Lines 147: computes `inter.zeeman.matrix` using `inter.zeeman.matrix=cell(1,nspins)`.
- Lines 148: computes `inter.zeeman.matrix{nspins}` using `inter.zeeman.matrix{nspins}=props.g_tensor.matrix`.
- Lines 153: computes `inter.coupling.matrix` using `inter.coupling.matrix=mat2cell(zeros(3*nspins,3*nspins),3*ones(nspins,1),3*ones(nspins,1))`.
- Lines 155: computes `inter.coupling.matrix{n,end}` using `inter.coupling.matrix{n,end}=1e6*gauss2mhz(props.hfc.full.matrix{index(n)}/2)`.
- Lines 156: computes `inter.coupling.matrix{end,n}` using `inter.coupling.matrix{end,n}=1e6*gauss2mhz(props.hfc.full.matrix{index(n)}/2)`.
- Lines 164: computes `inter.coupling.matrix{n,k}` using `inter.coupling.matrix{n,k}=[]`.
- Lines 172: computes `killing_pattern` using `killing_pattern=cellfun(@isempty,inter.coupling.matrix(:,nspins)); killing_pattern(end)=0`.
- Lines 173: computes `inter.coupling.matrix(:,killing_pattern)` using `inter.coupling.matrix(:,killing_pattern)=[]`.
- Lines 174: computes `inter.coupling.matrix(killing_pattern,:)` using `inter.coupling.matrix(killing_pattern,:)=[]`.
- Lines 175: computes `inter.zeeman.matrix(killing_pattern)` using `inter.zeeman.matrix(killing_pattern)=[]`.
- Lines 176: computes `sys.isotopes(killing_pattern)` using `sys.isotopes(killing_pattern)=[]`.
- Lines 190: computes `inter.zeeman.matrix{n}` using `inter.zeeman.matrix{n}=-props.cst{index(n)}+eye(3)*references(ref_index(n))`.

### Local helper functions

- Line 269: `grumble()` — `function grumble(props,nuclei,references,options)`.
  - Representative operation: `if ~isstruct(props)`.
  - Representative operation: `error('the first argument must be a structure returned by gparse().')`.

## Parameters / inputs

- props -the output of gparse() function
- particles -a cell array of the following form:
- {{'H','1H'},{'N','15N'}...}
- giving the list of elements and isotopes that
- should be imported. If the isotope list contains
- an electron, e.g. {{'E','E'},{'H','1H'}...},
- then EPR mode is assumed -chemical shielding
- and scalar couplings are ignored, but g-tensor
- and hyperfine couplings are included.
- references -a vector of absolute shielding values for
- the reference substances that are to be
- placed at zero ppm chemical shift; you ne-
- ed to run separate electronic structure
- theory calculations for those substances
- wuth the same method. Absolute isotropic
- shielding values for tetramethylsilane in
- vacuum are:
- GIAO 13C 1H
- B3LYP/6-31G* 189.6621 32.1833
- B3LYP/6-311+G(2d,p) 182.4485 31.8201
- HF/6-31G* 199.9711 32.5957
- HF/6-311+G(2d,p) 192.5828 32.0710
- CSGT 13C 1H
- B3LYP/6-31G* 188.5603 29.1952
- B3LYP/6-311+G(2d,p) 182.1386 31.7788
- HF/6-31G* 196.8670 29.5517
- HF/6-311+G(2d,p) 192.5701 31.5989
- This setting is ignored when electrons are
- present in the system.
- options.min_j -scalar coupling threshold in Hz. J-coup-
- lings smaller than this value will be
- ignored in the NMR mode.
- options.min_hfc -hyperfine coupling threshold in Hz. Hy-
- perfine tensors with a Frobenius norm
- smaller than this value will be ignored
- in the EPR mode.
- options.purge -if set to 'on' in EPR mode, removes the
- spins with hyperfine coupling below
- options.min_hfc from the spin system.
- options.no_xyz -if set to 1, causes the function to ig-
- nore the coordinate information and
- only keep the interaction tensors

## Outputs

- sys.isotopes Nspins x 1 cell array of strings
- inter.coordinates Nspins x 3 dense matrix, Angstrom. Not
- returned if there is an electron in the
- isotope list (in the EPR case it is not
- a good idea to use the molecular
- coordinates for spins).
- inter.zeeman.matrix Nspins x 1 cell array of 3x3 matrices,
- ppm for nuclei, g-tensor for electrons.
- Zero interactions have zero matrices.
- inter.coupling.matrix Nspins x Nspins cell array of 3x3 mat-
- rices, all in Hz. Zero interactions have
- zero matrices.
- inter.coupling.scalar Nspins x Nspins cell array of scalar
- couplings, all in Hz. Zero couplings are
- returned as zeros.
- inter.spinrot.matrix spin-rotation coupling tensors for
- each nucleus

## Implementation structure

- Makes Spinach data structures from parsed outputs of electronic
- structure theory packages, such as Gaussian and ORCA. Syntax:
- [sys,inter]=g2spinach(props,particles,references,options)
- props -the output of gparse() function
- particles -a cell array of the following form:
- {{'H','1H'},{'N','15N'}...}
- giving the list of elements and isotopes that
- should be imported. If the isotope list contains
- an electron, e.g. {{'E','E'},{'H','1H'}...},
- then EPR mode is assumed -chemical shielding
- and scalar couplings are ignored, but g-tensor
- and hyperfine couplings are included.

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `exist()`, `grumble()`, `strcmp()`, `true()`, `isfield()`, `false()`, `num2cell()`, `ismember()`, `mat2cell()`, `gauss2mhz()`, `index()`, `cellfun()`, `killing_pattern()`, `any()`, `references()`, `ref_index()`.
