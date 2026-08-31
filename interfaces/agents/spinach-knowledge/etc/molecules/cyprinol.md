# etc/molecules/cyprinol.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/etc/molecules/cyprinol.m`
- Signature: `[sys,inter,bas]=cyprinol()`
- Total lines: 193

## Purpose

Spin system of cyprinol. Isotropic chemical shifts and J-couplings are taken from http://dx.doi.org/10.1002/mrc.4782 and, when not gi- ven there, estimated by tossing a twenty-sided coin. Syntax: [sys,inter,bas]=cyprinol()

## Physical / mathematical content

- This file belongs to the `etc` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content


## Code-derived implementation details

### Comment-guided execution stages

- Lines 25-26: Hydrogen atoms; implemented by `H_iso=cell(1,42); H_iso(:)={'1H'}; numH=numel(H_iso)`.
- Lines 28-29: Shorthands for human-readable coupling designations below; implemented by `H1a=1; H1b=2; H2=2; H3=3; H4a=5; H4b=6; H5=7; H6a=8`.
- Lines 36-37: Proton chemical shifts; implemented by `H_CS=cell(1,numH)`.
- Lines 50-55: TODO: proton J-couplings; implemented by `C_iso=cell(1,27); C_iso(:)={'13C'}; numC=numel(C_iso)`.
- Lines 52-55: TODO: proton coordinates; implemented by `C_iso=cell(1,27); C_iso(:)={'13C'}; numC=numel(C_iso)`.
- Lines 54-55: Carbon atoms; implemented by `C_iso=cell(1,27); C_iso(:)={'13C'}; numC=numel(C_iso)`.
- Lines 57-58: Carbon chemical shifts; implemented by `C_CS=cell(1,numC)`.
- Lines 67-68: Carbon J-coupling estimates; implemented by `CC_J=cell(numC)`.
- Lines 109-112: TODO: carbon coordinates; implemented by `CH_J=cell(numC,numH)`.
- Lines 111-112: Carbon-proton J-coupling estimates; implemented by `CH_J=cell(numC,numH)`.
- Lines 169-170: Combine isotope arrays; implemented by `sys.isotopes=[H_iso C_iso]`.
- Lines 172-173: Combine chemical shift arrays; implemented by `inter.zeeman.scalar=[H_CS C_CS]`.
- Lines 175-176: Combine J-coupling arrays; implemented by `inter.coupling.scalar=cell(numH+numC)`.
- Lines 180-181: Symmetry settings; implemented by `bas.sym_group={'S3','S3','S3'}`.

### Key state/data transformations

- Lines 26: computes `H_iso` using `H_iso=cell(1,42); H_iso(:)={'1H'}; numH=numel(H_iso)`.
- Lines 29: computes `H1a` using `H1a=1; H1b=2; H2=2; H3=3; H4a=5; H4b=6; H5=7; H6a=8`.
- Lines 30: computes `H6b` using `H6b=9; H7=10; H8=11; H9=12; H11a=13; H11b=14; H12=15; H14=16`.
- Lines 31: computes `H15a` using `H15a=17; H15b=18; H16a=19; H16b=20; H17=21; H18a=22; H19a=23; H20=24`.
- Lines 32: computes `H21a` using `H21a=25; H22a=26; H22b=27; H23a=28; H23b=29; H24a=30; H24b=31; H25=32`.
- Lines 33: computes `H26a` using `H26a=33; H26b=34; H27a=35; H27b=36; H18b=37; H18c=38; H19b=39; H19c=40`.
- Lines 34: computes `H21b` using `H21b=41; H21c=42`.
- Lines 37: computes `H_CS` using `H_CS=cell(1,numH)`.
- Lines 38: computes `H_CS{H1a}` using `H_CS{H1a}= 1.42; H_CS{H17}= 1.85; H_CS{H1b}= 1.42; H_CS{H18a}=0.72`.
- Lines 39: computes `H_CS{H2}` using `H_CS{H2}= 1.64; H_CS{H19a}=0.82; H_CS{H3}= 3.99; H_CS{H20}= 1.40`.
- Lines 40: computes `H_CS{H4a}` using `H_CS{H4a}= 1.51; H_CS{H21a}=1.02; H_CS{H4b}= 1.32; H_CS{H22a}=1.46`.
- Lines 41: computes `H_CS{H5}` using `H_CS{H5}= 2.15; H_CS{H22b}=1.11; H_CS{H6a}= 1.42; H_CS{H23a}=1.48`.
- Lines 42: computes `H_CS{H6b}` using `H_CS{H6b}= 1.35; H_CS{H23b}=1.31; H_CS{H7}= 3.79; H_CS{H24a}=1.35`.
- Lines 43: computes `H_CS{H8}` using `H_CS{H8}= 1.48; H_CS{H24b}=1.35; H_CS{H9}= 1.68; H_CS{H25}= 1.80`.
- Lines 44: computes `H_CS{H11a}` using `H_CS{H11a}=1.68; H_CS{H26a}=4.05; H_CS{H11b}=1.57; H_CS{H26b}=4.00`.
- Lines 45: computes `H_CS{H12}` using `H_CS{H12}= 3.95; H_CS{H27a}=3.60; H_CS{H14}= 1.94; H_CS{H27b}=3.56`.
- Lines 46: computes `H_CS{H15a}` using `H_CS{H15a}=1.77; H_CS{H18b}=0.72; H_CS{H15b}=1.12; H_CS{H18c}=0.72`.
- Lines 47: computes `H_CS{H16a}` using `H_CS{H16a}=1.89; H_CS{H19b}=0.82; H_CS{H16b}=1.28; H_CS{H19c}=0.82`.

## Outputs

- sys -Spinach spin system description structure
- inter -Spinach interaction description structure
- bas -Spinach basis set description structure
- Note: if you are looking for a test spin system, strychnine.m is a
- more complete alternative.
- Bud MacAulay
- Ilya Kuprov

## Implementation structure

- Spin system of cyprinol. Isotropic chemical shifts and J-couplings
- are taken from http://dx.doi.org/10.1002/mrc.4782 and, when not gi-
- ven there, estimated by tossing a twenty-sided coin. Syntax:
- [sys,inter,bas]=cyprinol()
- sys -Spinach spin system description structure
- inter -Spinach interaction description structure
- bas -Spinach basis set description structure
- Note: if you are looking for a test spin system, strychnine.m is a
- more complete alternative.
- Bud MacAulay
- Ilya Kuprov
- Hydrogen atoms

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `H_iso()`, `C_iso()`.
