# examples/nmr_paramag/porphyrin_example_3.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_paramag/porphyrin_example_3.m`
- Signature: `porphyrin_example_3()`
- Total lines: 74

## Purpose

Computes PCS using different models in basic Cu(II) porphyrin complex. See the "getting started" manual at The paper describing the distributed PCS model used below is available at Calculation time: minutes, 64GB of RAM required

## Physical / mathematical content

- Paramagnetic NMR examples. These files work with pseudocontact shifts, paramagnetic relaxation, susceptibility tensors, and inverse problems for metal-site localisation or distributed electron density reconstruction.

## Numerical / algorithmic content


## Code-derived implementation details

### Comment-guided execution stages

- Lines 16-17: Porphyrin ring proton coordinates; implemented by `nxyz=[ 4.551635888 2.658552774 0.000000000`.
- Lines 30-31: Cu(II) g-tensor eigenvalues; implemented by `g_cu=diag([2.0000 2.0000 2.2000])`.
- Lines 33-34: Curie susceptibility tensor; implemented by `chi_cu=g2chi(g_cu,298,1/2)`.
- Lines 36-37: Parse ORCA log; implemented by `props=oparse('cu_porph_hfc.out')`.
- Lines 39-40: Extract hyperfine tensors; implemented by `hfcs=props.hfc.full.matrix(26:37)`.
- Lines 42-43: Compute HFC PCS; implemented by `hfc_pcs_cu=zeros(numel(hfcs),1)`.
- Lines 48-49: Parse ORCA cube and pad the density with zeros to avoid PBC effects; implemented by `pad_size=2; [density,ext]=ocparse('cu_porph_sd120.spindens.3d',pad_size)`.
- Lines 51-54: Compute extentens for original density without padding; implemented by `zoom_vol=[pad_size/(2*pad_size+1) (pad_size+1)/(2*pad_size+1) pad_size/(2*pad_size+1) (pad_size+1)/(2*pad_size+1) pad_size/(2*pad_size+1) (pad_size+1)/(2*pad_size+1)]`.
- Lines 56-57: Draw the probability density schematic; implemented by `[zoom_den,zoom_ext]=zoom_3d(density,ext,zoom_vol)`.
- Lines 61-62: Solve Kuprov equation; implemented by `[pde_pcs_cu,pcs_cube]=kpcs(density,chi_cu,ext,nxyz,'fft')`.
- Lines 64-65: Compare HFC PCS with PDE PCS; implemented by `disp('Pseudocontact shifts [hfc, pde], ppm')`.
- Lines 68-69: Plot the PCS field schematic; implemented by `[zoom_pcs,zoom_ext]=zoom_3d(pcs_cube,ext,zoom_vol)`.

### Control flow inferred from the code

- Line 44: `for` loop over `n=1:numel(hfcs)`.

### Key state/data transformations

- Lines 17: computes `nxyz` using `nxyz=[ 4.551635888 2.658552774 0.000000000`.
- Lines 31: computes `g_cu` using `g_cu=diag([2.0000 2.0000 2.2000])`.
- Lines 34: computes `chi_cu` using `chi_cu=g2chi(g_cu,298,1/2)`.
- Lines 37: computes `props` using `props=oparse('cu_porph_hfc.out')`.
- Lines 40: computes `hfcs` using `hfcs=props.hfc.full.matrix(26:37)`.
- Lines 43: computes `hfc_pcs_cu` using `hfc_pcs_cu=zeros(numel(hfcs),1)`.
- Lines 45: computes `hfc_pcs_cu(n)` using `hfc_pcs_cu(n)=hfc2pcs(hfcs{n},chi_cu,'1H')`.
- Lines 49: computes `pad_size` using `pad_size=2; [density,ext]=ocparse('cu_porph_sd120.spindens.3d',pad_size)`.
- Lines 52-54: computes `zoom_vol` using `zoom_vol=[pad_size/(2*pad_size+1) (pad_size+1)/(2*pad_size+1) pad_size/(2*pad_size+1) (pad_size+1)/(2*pad_size+1) pad_size/(2*pad_size+1) (pad_size+1)/(2*pad_size+1)]`.
- Lines 57: computes `[zoom_den,zoom_ext]` using `[zoom_den,zoom_ext]=zoom_3d(density,ext,zoom_vol)`.
- Lines 62: computes `[pde_pcs_cu,pcs_cube]` using `[pde_pcs_cu,pcs_cube]=kpcs(density,chi_cu,ext,nxyz,'fft')`.
- Lines 69: computes `[zoom_pcs,zoom_ext]` using `[zoom_pcs,zoom_ext]=zoom_3d(pcs_cube,ext,zoom_vol)`.

## Implementation structure

- Computes PCS using different models in basic Cu(II) porphyrin complex. See
- the "getting started" manual at
- The paper describing the distributed PCS model used below is available at
- Calculation time: minutes, 64GB of RAM required
- Porphyrin ring proton coordinates
- Cu(II) g-tensor eigenvalues
- Curie susceptibility tensor
- Parse ORCA log
- Extract hyperfine tensors
- Compute HFC PCS
- Parse ORCA cube and pad the density with zeros to avoid PBC effects
- Compute extentens for original density without padding

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `g2chi()`, `oparse()`, `hfc_pcs_cu()`, `hfc2pcs()`, `ocparse()`, `zoom_3d()`, `kfigure()`, `volplot()`, `molplot()`, `kpcs()`.
