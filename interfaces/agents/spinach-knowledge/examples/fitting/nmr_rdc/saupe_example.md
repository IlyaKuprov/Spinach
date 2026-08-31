# examples/fitting/nmr_rdc/saupe_example.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/fitting/nmr_rdc/saupe_example.m`
- Signature: `saupe_example()`
- Total lines: 50

## Purpose

Extracting Saupe order matrix from NH RDC data. Experimental measurements kindly provided by Andras Boeszoermenyi, Thibault Viennet, and Hari Arthanari.

## Physical / mathematical content

- Fitting examples. These files formulate parameter-estimation workflows in which simulated spectra or observables are matched to data, usually through nonlinear optimisation, residual construction, and physically constrained parameterisations.

## Numerical / algorithmic content


## Code-derived implementation details

### Comment-guided execution stages

- Lines 9-10: Read the PDB file; implemented by `[pdb_aa,~,pdb_id,coords]=read_pdb_pro('protein.pdb',1)`.
- Lines 12-14: Read RDC data; implemented by `load('rdc_data.mat','aa_num','rdc_hz', 'isotope_a','isotope_b','atom_a','atom_b')`.
- Lines 16-17: Make isotope table; implemented by `isotopes=[isotope_a isotope_b]`.
- Lines 19-20: Match up RDCs with coordinates; implemented by `for n=1:numel(rdc_hz)`.
- Lines 22-23: Locate both atoms; implemented by `index_a=(pdb_aa==aa_num(n))&strcmp(pdb_id,atom_a{n})`.
- Lines 26-27: Extract coordinates; implemented by `xyz{n,1}=coords{index_a}; xyz{n,2}=coords{index_b}`.
- Lines 31-32: Call RDC fitter; implemented by `S=rdc_fit(isotopes,xyz,rdc_hz)`.
- Lines 35-36: Back-calculate RDCs; implemented by `for n=1:numel(rdc_hz)`.
- Lines 41-42: Do the plotting; implemented by `kfigure(); plot(rdc_hz,rdc_theo','r.')`.

### Control flow inferred from the code

- Line 20: `for` loop over `n=1:numel(rdc_hz)`.
- Line 36: `for` loop over `n=1:numel(rdc_hz)`.

### Key state/data transformations

- Lines 10: computes `[pdb_aa,~,pdb_id,coords]` using `[pdb_aa,~,pdb_id,coords]=read_pdb_pro('protein.pdb',1)`.
- Lines 17: computes `isotopes` using `isotopes=[isotope_a isotope_b]`.
- Lines 27: computes `xyz{n,1}` using `xyz{n,1}=coords{index_a}; xyz{n,2}=coords{index_b}`.
- Lines 32: computes `S` using `S=rdc_fit(isotopes,xyz,rdc_hz)`.
- Lines 37-38: computes `rdc_theo(n)` using `rdc_theo(n)=xyz2rdc(isotope_a{n},isotope_b{n}, xyz{n,1},xyz{n,2},{S,'saupe'})`.

## Implementation structure

- Extracting Saupe order matrix from NH RDC data. Experimental
- measurements kindly provided by Andras Boeszoermenyi, Thibault
- Viennet, and Hari Arthanari.
- Read the PDB file
- Read RDC data
- Make isotope table
- Match up RDCs with coordinates
- Locate both atoms
- Extract coordinates
- Call RDC fitter
- Back-calculate RDCs
- Do the plotting

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `read_pdb_pro()`, `load()`, `aa_num()`, `strcmp()`, `rdc_fit()`, `rdc_theo()`, `xyz2rdc()`, `kfigure()`, `kxlabel()`, `kylabel()`, `xlim()`, `ylim()`.
