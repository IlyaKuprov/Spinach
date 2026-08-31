# examples/nmr_paramag/porphyrin_example_1.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_paramag/porphyrin_example_1.m`
- Signature: `porphyrin_example_1()`
- Total lines: 46

## Purpose

Computing PCS using different models in basic Cu(II) and Co(II) porphyrin complexes. The metal is at the origin. See the "getting started" manual at

## Physical / mathematical content

- Paramagnetic NMR examples. These files work with pseudocontact shifts, paramagnetic relaxation, susceptibility tensors, and inverse problems for metal-site localisation or distributed electron density reconstruction.

## Numerical / algorithmic content


## Code-derived implementation details

### Comment-guided execution stages

- Lines 10-11: Porphyrin ring proton coordinates; implemented by `nxyz=[ 4.551635888 2.658552774 0.000000000`.
- Lines 24-25: Co(II) g-tensor; implemented by `g_co=diag([3.0 3.0 2.0])`.
- Lines 27-28: Cu(II) g-tensor; implemented by `g_cu=diag([2.0 2.0 2.2])`.
- Lines 30-31: Curie susceptibility tensors; implemented by `chi_co=g2chi(g_co,298,1/2)`.
- Lines 34-35: Metal position; implemented by `mxyz=[0 0 0]`.
- Lines 37-38: PCS calculation; implemented by `point_pcs_co=ppcs(nxyz,mxyz,chi_co)`.
- Lines 41-42: Output; implemented by `disp('PCS [Co, Cu], ppm')`.

### Key state/data transformations

- Lines 11: computes `nxyz` using `nxyz=[ 4.551635888 2.658552774 0.000000000`.
- Lines 25: computes `g_co` using `g_co=diag([3.0 3.0 2.0])`.
- Lines 28: computes `g_cu` using `g_cu=diag([2.0 2.0 2.2])`.
- Lines 31: computes `chi_co` using `chi_co=g2chi(g_co,298,1/2)`.
- Lines 32: computes `chi_cu` using `chi_cu=g2chi(g_cu,298,1/2)`.
- Lines 35: computes `mxyz` using `mxyz=[0 0 0]`.
- Lines 38: computes `point_pcs_co` using `point_pcs_co=ppcs(nxyz,mxyz,chi_co)`.
- Lines 39: computes `point_pcs_cu` using `point_pcs_cu=ppcs(nxyz,mxyz,chi_cu)`.

## Implementation structure

- Computing PCS using different models in basic Cu(II) and Co(II) porphyrin
- complexes. The metal is at the origin. See the "getting started" manual at
- Porphyrin ring proton coordinates
- Co(II) g-tensor
- Cu(II) g-tensor
- Curie susceptibility tensors
- Metal position
- PCS calculation
- Output

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `g2chi()`, `ppcs()`.
