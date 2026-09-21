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
