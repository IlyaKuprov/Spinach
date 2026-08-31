# examples/nmr_paramag/dft_density.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_paramag/dft_density.m`
- Signature: `dft_density()`
- Total lines: 72

## Purpose

Simulation of the pseudocontact shift field of the Europium(III) complex of 1,4,7,10-tetrakis(2-pyridylmethyl)-1,4,7,10-tetraazacyclododecane. The spin density, the hyperfine couplings and the susceptibility tensor are imported from a DFT calculation. The partial differential equation used for the delo- calised model solution is described in: One outlier point is due to the presence of contact shift due to the isotro

## Physical / mathematical content

- Paramagnetic NMR examples. These files work with pseudocontact shifts, paramagnetic relaxation, susceptibility tensors, and inverse problems for metal-site localisation or distributed electron density reconstruction.

## Numerical / algorithmic content


## Code-derived implementation details

### Comment-guided execution stages

- Lines 18-19: Load unpaired electron probability density; implemented by `load('tetra_py_probden.mat','ext','xyz','probden','dx')`.
- Lines 21-22: Load DFT data (HFCs are read in Gauss); implemented by `props=gparse('tetra_py_dft_run.log')`.
- Lines 24-25: Normalize probability density; implemented by `probden=probden/(trapz(trapz(trapz(probden)))*dx^3)`.
- Lines 27-28: Get susceptibility tensor; implemented by `[~,~,rank2]=mat2sphten(props.chi)`.
- Lines 31-32: Get point pseudocontact shifts; implemented by `pcs_point=ppcs(xyz,[0 0 0],chi)`.
- Lines 34-35: Solve Kuprov equation; implemented by `[pcs_distr,pcs_3d]=kpcs(probden,chi,ext,xyz,'fft')`.
- Lines 37-38: Get DFT pseudocontact shifts; implemented by `for n=1:(props.natoms-1)`.
- Lines 40-41: Assign isotopes; implemented by `if strcmp(props.symbols{n},'H')`.
- Lines 49-50: Compute DFT PCS; implemented by `pcs_hfc(n)=hfc2pcs(props.hfc.full.matrix{n},chi,isotope)`.
- Lines 54-55: Plot DFT against distributed; implemented by `kfigure(); plot(pcs_hfc,pcs_distr,'ro')`.
- Lines 60-61: Plot point against distributed; implemented by `kfigure(); plot(pcs_point,pcs_distr,'ro')`.
- Lines 66-67: Plot the distributed solution; implemented by `kfigure(); pcs_3d=sqrt(abs(pcs_3d)).*sign(pcs_3d)`.

### Control flow inferred from the code

- Line 38: `for` loop over `n=1:(props.natoms-1)`.
- Line 41: conditional branch on `strcmp(props.symbols{n},'H')`.

### Key state/data transformations

- Lines 22: computes `props` using `props=gparse('tetra_py_dft_run.log')`.
- Lines 25: computes `probden` using `probden=probden/(trapz(trapz(trapz(probden)))*dx^3)`.
- Lines 28: computes `[~,~,rank2]` using `[~,~,rank2]=mat2sphten(props.chi)`.
- Lines 29: computes `chi` using `chi=sphten2mat(0,[0 0 0],rank2)`.
- Lines 32: computes `pcs_point` using `pcs_point=ppcs(xyz,[0 0 0],chi)`.
- Lines 35: computes `[pcs_distr,pcs_3d]` using `[pcs_distr,pcs_3d]=kpcs(probden,chi,ext,xyz,'fft')`.
- Lines 42: computes `isotope` using `isotope='1H'`.
- Lines 50: computes `pcs_hfc(n)` using `pcs_hfc(n)=hfc2pcs(props.hfc.full.matrix{n},chi,isotope)`.
- Lines 67: computes `kfigure(); pcs_3d` using `kfigure(); pcs_3d=sqrt(abs(pcs_3d)).*sign(pcs_3d)`.

## Implementation structure

- Simulation of the pseudocontact shift field of the Europium(III) complex of
- 1,4,7,10-tetrakis(2-pyridylmethyl)-1,4,7,10-tetraazacyclododecane. The spin
- density, the hyperfine couplings and the susceptibility tensor are imported
- from a DFT calculation. The partial differential equation used for the delo-
- calised model solution is described in:
- One outlier point is due to the presence of contact shift due to the isotro-
- pic hyperfine coupling being non-sero for that particular nucleus. Point mo-
- del and Kuprov equation do not include contact shifts.
- Load unpaired electron probability density
- Load DFT data (HFCs are read in Gauss)
- Normalize probability density
- Get susceptibility tensor

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `load()`, `gparse()`, `trapz()`, `mat2sphten()`, `sphten2mat()`, `ppcs()`, `kpcs()`, `strcmp()`, `pcs_hfc()`, `hfc2pcs()`, `kfigure()`, `kxlabel()`, `kylabel()`, `sign()`, `volplot()`, `molplot()`.
