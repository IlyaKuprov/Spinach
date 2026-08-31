# etc/diamond_defects/diamond_p1_13c.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/etc/diamond_defects/diamond_p1_13c.m`
- Signature: `[sys,inter]=diamond_p1_13c(parameters)`
- Total lines: 330

## Purpose

P1 centre spin system with 13C neighbours in diamond. Syntax: [sys,inter]=diamond_p1_13c(parameters) The electron and nitrogen parameters follow diamond_p1.m. The 13C hyperfine tensors and site assignments are from: R C Barklie and J Guven, J. Phys. C: Solid State Phys. 14, 3621-3631 (1981), doi:10.1088/0022-3719/14/25/009 A Cox, M E Newton, and J M Baker, J. Phys.: Condens. Matter 6, 551-563 (1994), doi:10.1088/0953

## Physical / mathematical content

- This file belongs to the `etc` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.
- The optimisation logic is Newton or Newton-like: search directions use first- and second-order local curvature information, usually with regularisation or line-search safeguards.
- The spin physics includes through-space magnetic dipole-dipole coupling, a rank-2 anisotropic interaction with strong orientation dependence and characteristic secular/non-secular structure.
- Quadrupolar physics is relevant: nuclei with spin > 1/2 interact with the electric field gradient tensor, introducing second-rank anisotropy, asymmetry, and overtone or MQ phenomena.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 53-54: Check input count; implemented by `if nargin~=1`.
- Lines 58-59: Check consistency; implemented by `grumble(parameters)`.
- Lines 61-62: Set crystal-to-defect frame rotation; implemented by `crys2def=rotmat_align([1 1 1],[0 0 1])`.
- Lines 64-65: Set the same orientation convention as diamond_p1.m; implemented by `switch parameters.orientation`.
- Lines 81-82: Complain and bomb out; implemented by `error('unknown orientation specification.')`.
- Lines 86-90: Define the 13C site labels; implemented by `site_lbl={'G1_C1','G2_C3','G3_calc','G4_C5','G5_calc', 'G6_calc','G7_calc','G8_C2','G9_calc','G10_calc', 'G11_calc','G12_calc','G13_calc','G14_C4', 'G15_calc','G16_calc'…`.
- Lines 92-93: Define spin system isotopes and labels; implemented by `sys.isotopes=[{'E',parameters.nitrogen},repmat({'13C'},1,numel(site_lbl))]`.
- Lines 98-99: Electron g-tensor; implemented by `inter.zeeman.matrix{1}=R*diag([2.00220 2.00220 2.00218])*R'`.
- Lines 101-102: HFC and NQI tensor; implemented by `switch parameters.nitrogen`.
- Lines 106-107: Hyperfine and quadrupolar coupling; implemented by `inter.coupling.matrix{1,2}=R*diag([81.3e6 81.3e6 114.0e6])*R'`.
- Lines 112-113: Only hyperfine coupling, opposite sign; implemented by `inter.coupling.matrix{1,2}=R*diag([-114.0e6 -114.0e6 -159.9e6])*R'`.
- Lines 117-118: Complain and bomb out; implemented by `error('wrong nitrogen isotope.')`.
- Lines 122-144: 13C site source map, in row order below: G1_C1 -C1, Cox et al. 1994, Table 2 G2_C3 -C3, Cox et al. 1994, Table 2, with sign from Peaker et al. 2016 G3_calc -Peaker et al. 2016, Table 4 G4_C5 -C5, Cox et al. 1994, Table 2, assignment from Peaker et al. 2016 G5_calc -Peaker et al. 2016, Table 4 G6_calc -Peaker et al. 2016, Table 4 G7_calc -Peaker et al. 2016, Table 4 G8_C2 -C2, Cox et al. 1994, Table 2, assignment from Peaker et al. 2016 G9_calc -Peaker et al. 2016, Table 4 G10_calc -Peaker et al. 2016, Table 4 G11_calc -Peaker et al. 2016, Table 4 G12_calc -Peaker et al. 2016, Table 4 G13_calc -Peaker et al. 2016, Table 4 G14_C4 -C4, Cox et al. 1994, Table 2, assignment from Peaker et al. 2016 G15_calc -Peaker et al. 2016, Table 4 G16_calc -Peaker et al. 2016, Table 4 G17_calc -Peaker et al. 2016, Table 4 G18_calc -Peaker et al. 2016, Table 4; implemented by `hfc_vals=[+139.531 +139.531 +338.171`.
- Lines 143-144: 13C hyperfine eigenvalues in MHz; implemented by `hfc_vals=[+139.531 +139.531 +338.171`.
- Lines 163-164: 13C hyperfine principal-axis polar angles, degrees; implemented by `hfc_theta=[90.00 35.264 54.736`.
- Lines 183-184: 13C hyperfine principal-axis azimuths, degrees; implemented by `hfc_phi=[135.00 225.00 45.00`.
- Lines 203-204: Representative ideal-lattice directions for G1...G18; implemented by `site_frac=[+0.25 +0.25 +0.25`.
- Lines 223-225: Peaker et al. tabulated distances from N-G1 bond mid-point, units of a0; implemented by `site_rad=[0.284 0.545 0.553 0.735 0.735 0.894 0.894 0.899 0.906, 1.024 1.033 1.242 1.242 1.246 1.252 1.514 1.596 1.947]`.

### Control flow inferred from the code

- Line 54: conditional branch on `nargin~=1`.
- Line 65: dispatches on `parameters.orientation`; cases `'100'`, `'110'`, `'111'`.
- Line 102: dispatches on `parameters.nitrogen`; cases `'14N'`, `'15N'`.
- Line 235: `for` loop over `n=1:numel(site_lbl)`.
- Line 258: `for` loop over `n=1:numel(site_lbl)`.
- Line 263: `for` loop over `n=1:numel(site_lbl)`.

### Key state/data transformations

- Lines 62: computes `crys2def` using `crys2def=rotmat_align([1 1 1],[0 0 1])`.
- Lines 69: computes `R` using `R=rotmat_align([1 1 1],[1 0 0])`.
- Lines 87-90: computes `site_lbl` using `site_lbl={'G1_C1','G2_C3','G3_calc','G4_C5','G5_calc', 'G6_calc','G7_calc','G8_C2','G9_calc','G10_calc', 'G11_calc','G12_calc','G13_calc','G14_C4', 'G15_calc','G16_calc'…`.
- Lines 93: computes `sys.isotopes` using `sys.isotopes=[{'E',parameters.nitrogen},repmat({'13C'},1,numel(site_lbl))]`.
- Lines 94: computes `sys.labels` using `sys.labels=[{'E',['N_' parameters.nitrogen]},site_lbl]`.
- Lines 95: computes `inter.zeeman.matrix` using `inter.zeeman.matrix=cell(1,numel(sys.isotopes))`.
- Lines 96: computes `inter.coupling.matrix` using `inter.coupling.matrix=cell(numel(sys.isotopes),numel(sys.isotopes))`.
- Lines 99: computes `inter.zeeman.matrix{1}` using `inter.zeeman.matrix{1}=R*diag([2.00220 2.00220 2.00218])*R'`.
- Lines 107: computes `inter.coupling.matrix{1,2}` using `inter.coupling.matrix{1,2}=R*diag([81.3e6 81.3e6 114.0e6])*R'`.
- Lines 108: computes `inter.coupling.matrix{2,2}` using `inter.coupling.matrix{2,2}=R*zfs2mat(-3.97e6,0,0,0,0)*R'`.
- Lines 144: computes `hfc_vals` using `hfc_vals=[+139.531 +139.531 +338.171`.
- Lines 164: computes `hfc_theta` using `hfc_theta=[90.00 35.264 54.736`.
- Lines 184: computes `hfc_phi` using `hfc_phi=[135.00 225.00 45.00`.
- Lines 204: computes `site_frac` using `site_frac=[+0.25 +0.25 +0.25`.
- Lines 224-225: computes `site_rad` using `site_rad=[0.284 0.545 0.553 0.735 0.735 0.894 0.894 0.899 0.906, 1.024 1.033 1.242 1.242 1.246 1.252 1.514 1.596 1.947]`.
- Lines 228: computes `a0` using `a0=3.567`.
- Lines 229: computes `mid_frac` using `mid_frac=0.284*[1 1 1]/sqrt(3)`.
- Lines 232: computes `xyz_cub` using `xyz_cub=zeros(numel(site_lbl),3)`.

### Local helper functions

- Line 283: `grumble()` — `function grumble(parameters)`.
  - Representative operation: `if(~isstruct(parameters))`.
  - Representative operation: `error('parameters must be a structure.')`.

## Parameters / inputs

- a structure (parameters.*) with the following fields:
- .orientation -'111', '110', or '100' crystal plane normal
- aligned with the magnetic field
- .nitrogen -'14N' or '15N'

## Outputs

- sys -Spinach system specification structure
- inter -Spinach interaction specification structure

## Implementation structure

- P1 centre spin system with 13C neighbours in diamond. Syntax:
- [sys,inter]=diamond_p1_13c(parameters)
- The electron and nitrogen parameters follow diamond_p1.m. The 13C
- hyperfine tensors and site assignments are from:
- R C Barklie and J Guven, J. Phys. C: Solid State Phys. 14,
- 3621-3631 (1981), doi:10.1088/0022-3719/14/25/009
- A Cox, M E Newton, and J M Baker, J. Phys.: Condens. Matter
- 6, 551-563 (1994), doi:10.1088/0953-8984/6/2/012
- C V Peaker, M K Atumi, J P Goss, P R Briddon, A B Horsfall,
- M J Rayson, and R Jones, Diamond Relat. Mater. 70, 118-123
- (2016), doi:10.1016/j.diamond.2016.10.013
- Coordinates are representative site coordinates reconstructed from

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `rotmat_align()`, `zfs2mat()`, `site_frac()`, `xyz_cub()`, `site_rad()`, `xyz_lab()`, `sind()`, `hfc_theta()`, `cosd()`, `hfc_phi()`, `hfc_vals()`, `isstruct()`, `isfield()`, `ischar()`, `ismember()`.
