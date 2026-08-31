# etc/molecules/dac_reaction.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/etc/molecules/dac_reaction.m`
- Signature: `[sys,inter,bas,kin]=dac_reaction()`
- Total lines: 194

## Purpose

Example Diels-Alder cycloaddition reaction settings: pentadiene (reactant), acrylonitrile (reactant), exo-norbornene (product), endo-norbornene (product), and acetonitrile (solvent). Atom co- ordinates and chemical shift anisotropies are pulled from a DFT calculation; isotropic chemical shifts and J-couplings are ex- perimental (some J-coupling signs are missing). Syntax: [sys,inter,bas,kin]=dac_reaction()

## Physical / mathematical content

- This file belongs to the `etc` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.
- The relaxation model is Redfield-type perturbation theory: fluctuating interactions enter through correlation functions or spectral densities and generate a linear relaxation superoperator.

## Numerical / algorithmic content


## Code-derived implementation details

### Comment-guided execution stages

- Lines 31-32: Find own location; implemented by `own_folder=mfilename('fullpath')`.
- Lines 35-36: Cyclopentadiene (substance A); implemented by `props_a=gparse([own_folder 'cyclopentadiene.log'])`.
- Lines 39-41: Replace isotropic shifts; implemented by `inter_a.zeeman.matrix=shift_iso(inter_a.zeeman.matrix,[1 2 3 4 5 6], [6.628 6.427 6.427 6.628 2.797 2.797])`.
- Lines 43-44: Create labels; implemented by `sys_a.labels={'HB','HA','HAp','HBp','HCp','HC'}`.
- Lines 46-47: Replace J-couplings; implemented by `inter_a.coupling.scalar=zeros(6,6)`.
- Lines 64-65: Acrylonitrile (substance B); implemented by `props_b=gparse([own_folder 'acrylonitrile.log'])`.
- Lines 68-70: Replace isotropic shifts; implemented by `inter_b.zeeman.matrix=shift_iso(inter_b.zeeman.matrix,[1 2 3], [6.203 6.077 5.703])`.
- Lines 72-73: Create labels; implemented by `sys_b.labels={'H7','H8','H9'}`.
- Lines 75-76: Replace couplings; implemented by `inter_b.coupling.scalar=zeros(3,3)`.
- Lines 82-83: Endo-norbornene carbonitrile (substance C); implemented by `props_c=gparse([own_folder 'norbornene_endo.log'])`.
- Lines 86-89: Replace isotropic shifts; implemented by `inter_c.zeeman.matrix=shift_iso(inter_c.zeeman.matrix,[1 2 3 4 5 6 7 8 9], [6.262 6.129 2.074 1.263 2.781 3.167 2.959 1.133 1.449])`.
- Lines 91-92: Create labels; implemented by `sys_c.labels={'H10','H11','H14','H15','H13','H12','H16','H17','H18'}`.
- Lines 94-95: Replace couplings (signs to be added); implemented by `inter_c.coupling.scalar=zeros(9,9)`.
- Lines 112-113: Exo-norbornene carbonitrile (substance D); implemented by `props_d=gparse([own_folder 'norbornene_exo.log'])`.
- Lines 116-119: Replace isotropic shifts; implemented by `inter_d.zeeman.matrix=shift_iso(inter_d.zeeman.matrix,[1 2 3 4 5 6 7 8 9], [6.099 5.9713 1.8972 1.4956 3.158 2.9830 1.4958 1.4958 2.113])`.
- Lines 121-122: Create labels; implemented by `sys_d.labels={'H19','H20','H23','H24','H21','H25','H26','H27','H22'}`.
- Lines 124-125: Replace couplings (signs to be added); implemented by `inter_d.coupling.scalar=zeros(9,9)`.
- Lines 142-143: Acetonitrile (substance E, solvent); implemented by `sys_e.isotopes={'1H','1H','1H'}`.

### Key state/data transformations

- Lines 32: computes `own_folder` using `own_folder=mfilename('fullpath')`.
- Lines 36: computes `props_a` using `props_a=gparse([own_folder 'cyclopentadiene.log'])`.
- Lines 37: computes `[sys_a,inter_a]` using `[sys_a,inter_a]=g2spinach(props_a,{{'H','1H'}},31.8)`.
- Lines 40-41: computes `inter_a.zeeman.matrix` using `inter_a.zeeman.matrix=shift_iso(inter_a.zeeman.matrix,[1 2 3 4 5 6], [6.628 6.427 6.427 6.628 2.797 2.797])`.
- Lines 44: computes `sys_a.labels` using `sys_a.labels={'HB','HA','HAp','HBp','HCp','HC'}`.
- Lines 47: computes `inter_a.coupling.scalar` using `inter_a.coupling.scalar=zeros(6,6)`.
- Lines 65: computes `props_b` using `props_b=gparse([own_folder 'acrylonitrile.log'])`.
- Lines 66: computes `[sys_b,inter_b]` using `[sys_b,inter_b]=g2spinach(props_b,{{'H','1H'}},31.8)`.
- Lines 69-70: computes `inter_b.zeeman.matrix` using `inter_b.zeeman.matrix=shift_iso(inter_b.zeeman.matrix,[1 2 3], [6.203 6.077 5.703])`.
- Lines 73: computes `sys_b.labels` using `sys_b.labels={'H7','H8','H9'}`.
- Lines 76: computes `inter_b.coupling.scalar` using `inter_b.coupling.scalar=zeros(3,3)`.
- Lines 77: computes `inter_b.coupling.scalar(idxof(sys_b,'H7'),idxof(sys_b,'H8'))` using `inter_b.coupling.scalar(idxof(sys_b,'H7'),idxof(sys_b,'H8'))= 0.91`.
- Lines 78: computes `inter_b.coupling.scalar(idxof(sys_b,'H7'),idxof(sys_b,'H9'))` using `inter_b.coupling.scalar(idxof(sys_b,'H7'),idxof(sys_b,'H9'))=17.92`.
- Lines 79: computes `inter_b.coupling.scalar(idxof(sys_b,'H8'),idxof(sys_b,'H9'))` using `inter_b.coupling.scalar(idxof(sys_b,'H8'),idxof(sys_b,'H9'))=11.75`.
- Lines 83: computes `props_c` using `props_c=gparse([own_folder 'norbornene_endo.log'])`.
- Lines 84: computes `[sys_c,inter_c]` using `[sys_c,inter_c]=g2spinach(props_c,{{'H','1H'}},31.8)`.
- Lines 87-89: computes `inter_c.zeeman.matrix` using `inter_c.zeeman.matrix=shift_iso(inter_c.zeeman.matrix,[1 2 3 4 5 6 7 8 9], [6.262 6.129 2.074 1.263 2.781 3.167 2.959 1.133 1.449])`.
- Lines 92: computes `sys_c.labels` using `sys_c.labels={'H10','H11','H14','H15','H13','H12','H16','H17','H18'}`.

## Parameters / inputs

- none

## Outputs

- sys, inter, bas -Spinach input data structures, remember
- to specify the field in sys.magnet
- kin -matching tables for which nuclei go where in which
- of the two chemical reactions

## Implementation structure

- Example Diels-Alder cycloaddition reaction settings: pentadiene
- (reactant), acrylonitrile (reactant), exo-norbornene (product),
- endo-norbornene (product), and acetonitrile (solvent). Atom co-
- ordinates and chemical shift anisotropies are pulled from a DFT
- calculation; isotropic chemical shifts and J-couplings are ex-
- perimental (some J-coupling signs are missing). Syntax:
- [sys,inter,bas,kin]=dac_reaction()
- none
- sys, inter, bas -Spinach input data structures, remember
- to specify the field in sys.magnet
- kin -matching tables for which nuclei go where in which
- of the two chemical reactions

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `mfilename()`, `own_folder()`, `gparse()`, `g2spinach()`, `shift_iso()`, `idxof()`, `num2cell()`, `merge_inp()`.
