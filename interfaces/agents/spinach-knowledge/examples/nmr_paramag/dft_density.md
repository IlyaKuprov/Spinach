# examples/nmr_paramag/dft_density.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_paramag/dft_density.m`
- Signature: `dft_density()`
- Total lines: 72

## Purpose

Simulation of the pseudocontact shift field of the Europium(III) complex of 1,4,7,10-tetrakis(2-pyridylmethyl)-1,4,7,10-tetraazacyclododecane. The spin density, the hyperfine couplings and the susceptibility tensor are imported from a DFT calculation. The partial differential equation used for the delo- calised model solution is described in: One outlier point is due to the presence of contact shift due to the isotro

## Physical / mathematical content

- Paramagnetic NMR examples. These files work with pseudocontact shifts, paramagnetic relaxation, susceptibility tensors, and inverse problems for metal-site localisation or distributed electron density reconstruction.

## Numerical / algorithmic content

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
