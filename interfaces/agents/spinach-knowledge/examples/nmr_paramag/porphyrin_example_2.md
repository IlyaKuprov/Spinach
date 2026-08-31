# examples/nmr_paramag/porphyrin_example_2.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_paramag/porphyrin_example_2.m`
- Signature: `porphyrin_example_2()`
- Total lines: 76

## Purpose

Computing PCS using different models in a basic Cu(II) porphyrin complex. The metal is at the origin. See the "getting started" manual at The paper describing the distributed PCS model used below is available at

## Physical / mathematical content

- Paramagnetic NMR examples. These files work with pseudocontact shifts, paramagnetic relaxation, susceptibility tensors, and inverse problems for metal-site localisation or distributed electron density reconstruction.

## Numerical / algorithmic content


## Code-derived implementation details

### Comment-guided execution stages

- Lines 14-15: Porphyrin ring proton coordinates; implemented by `nxyz=[ 4.551635888 2.658552774 0.000000000`.
- Lines 28-29: Cu(II) g-tensor; implemented by `g_cu=diag([2.0000 2.0000 2.2000])`.
- Lines 31-32: Curie susceptibility tensor; implemented by `chi_cu=g2chi(g_cu,298,1/2)`.
- Lines 34-35: Metal position; implemented by `mxyz=[0 0 0]`.
- Lines 37-38: PCS calculation using the point model; implemented by `point_pcs_cu=ppcs(nxyz,mxyz,chi_cu)`.
- Lines 40-41: Coordinates of all atoms with significant spin population; implemented by `xyz=[ 0.000000000000 0.000000000000 0.000000000000`.
- Lines 47-48: Mulliken spin populations; implemented by `rho=[0.6 0.1 0.1 0.1 0.1]'`.
- Lines 50-51: Multipole ranks to include; implemented by `L=0:14`.
- Lines 53-54: Multipole moments; implemented by `Ilm=points2mult(xyz,mxyz,rho,L,'points')`.
- Lines 56-57: PCS calculation using the distributed model; implemented by `distr_pcs_cu=lpcs(nxyz,mxyz,L,Ilm,chi_cu)`.
- Lines 59-60: Parse the ORCA log; implemented by `props=oparse('cu_porph_hfc.out')`.
- Lines 62-63: Extract hyperfine tensors; implemented by `hfcs=props.hfc.full.matrix(26:37)`.
- Lines 65-66: Compute DFT PCS; implemented by `dft_pcs_cu=zeros(size(point_pcs_cu))`.
- Lines 71-72: Comparison of point with distributed with DFT; implemented by `disp('PCS [point, distr, dft], ppm')`.

### Control flow inferred from the code

- Line 67: `for` loop over `n=1:numel(hfcs)`.

### Key state/data transformations

- Lines 15: computes `nxyz` using `nxyz=[ 4.551635888 2.658552774 0.000000000`.
- Lines 29: computes `g_cu` using `g_cu=diag([2.0000 2.0000 2.2000])`.
- Lines 32: computes `chi_cu` using `chi_cu=g2chi(g_cu,298,1/2)`.
- Lines 35: computes `mxyz` using `mxyz=[0 0 0]`.
- Lines 38: computes `point_pcs_cu` using `point_pcs_cu=ppcs(nxyz,mxyz,chi_cu)`.
- Lines 41: computes `xyz` using `xyz=[ 0.000000000000 0.000000000000 0.000000000000`.
- Lines 48: computes `rho` using `rho=[0.6 0.1 0.1 0.1 0.1]'`.
- Lines 51: computes `L` using `L=0:14`.
- Lines 54: computes `Ilm` using `Ilm=points2mult(xyz,mxyz,rho,L,'points')`.
- Lines 57: computes `distr_pcs_cu` using `distr_pcs_cu=lpcs(nxyz,mxyz,L,Ilm,chi_cu)`.
- Lines 60: computes `props` using `props=oparse('cu_porph_hfc.out')`.
- Lines 63: computes `hfcs` using `hfcs=props.hfc.full.matrix(26:37)`.
- Lines 66: computes `dft_pcs_cu` using `dft_pcs_cu=zeros(size(point_pcs_cu))`.
- Lines 68: computes `dft_pcs_cu(n)` using `dft_pcs_cu(n)=hfc2pcs(hfcs{n},chi_cu,'1H')`.

## Implementation structure

- Computing PCS using different models in a basic Cu(II) porphyrin complex.
- The metal is at the origin. See the "getting started" manual at
- The paper describing the distributed PCS model used below is available at
- Porphyrin ring proton coordinates
- Cu(II) g-tensor
- Curie susceptibility tensor
- Metal position
- PCS calculation using the point model
- Coordinates of all atoms with significant spin population
- Mulliken spin populations
- Multipole ranks to include
- Multipole moments

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `g2chi()`, `ppcs()`, `points2mult()`, `lpcs()`, `oparse()`, `dft_pcs_cu()`, `hfc2pcs()`.
