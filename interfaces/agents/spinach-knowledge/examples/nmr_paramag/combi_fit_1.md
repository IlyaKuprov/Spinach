# examples/nmr_paramag/combi_fit_1.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_paramag/combi_fit_1.m`
- Signature: `combi_fit_1()`
- Total lines: 50

## Purpose

Extracting the susceptibility tensor from DFT hyperfine tensors and experimental paramagnetic shifts. A combinatorial procedure is used that cycles through ambiguous assignments. Calculation time: minutes.

## Physical / mathematical content

- Paramagnetic NMR examples. These files work with pseudocontact shifts, paramagnetic relaxation, susceptibility tensors, and inverse problems for metal-site localisation or distributed electron density reconstruction.

## Numerical / algorithmic content


## Code-derived implementation details

### Comment-guided execution stages

- Lines 12-13: Read DFT HFCs in Gauss; implemented by `props=gparse('l2_parker_funk_tm.log')`.
- Lines 16-17: Isotope list; implemented by `parameters.isotopes=cell(27,1)`.
- Lines 22-23: Spin groups with identical PCS; implemented by `parameters.spin_groups={[28 22 30]; [25 33 81]; [80 32 24]; [27 82 21]`.
- Lines 26-27: Diamagnetic shifts; implemented by `parameters.d_shifts=[3.62 2.65 2.65 2.86 4.95 4.10 7.40 8.00 7.80]`.
- Lines 29-30: Diamagnetic shift ambiguities; implemented by `parameters.d_ambig={[1 2 3 4],[5 6],[7 8 9]}`.
- Lines 32-33: Paramagnetic shifts; implemented by `parameters.p_shifts=[6.5 -56.9 -21.7 -23.0 11.4 54.6 18.6 16.4 17.8]`.
- Lines 35-36: Paramagnetic shift ambiguities; implemented by `parameters.p_ambig={[3 4],[7 8 9]}`.
- Lines 38-39: Run the combinatorial fitting; implemented by `[~,~,pcs_theo,pcs_expt,~,~]=pcs_combi_fit(parameters)`.
- Lines 41-42: Plot theory against experiment; implemented by `kfigure()`.

### Control flow inferred from the code

- Line 18: `for` loop over `n=1:27`.

### Key state/data transformations

- Lines 13: computes `props` using `props=gparse('l2_parker_funk_tm.log')`.
- Lines 14: computes `parameters.hfcs` using `parameters.hfcs=props.hfc.full.matrix`.
- Lines 17: computes `parameters.isotopes` using `parameters.isotopes=cell(27,1)`.
- Lines 19: computes `parameters.isotopes{n}` using `parameters.isotopes{n}='1H'`.
- Lines 23: computes `parameters.spin_groups` using `parameters.spin_groups={[28 22 30]; [25 33 81]; [80 32 24]; [27 82 21]`.
- Lines 27: computes `parameters.d_shifts` using `parameters.d_shifts=[3.62 2.65 2.65 2.86 4.95 4.10 7.40 8.00 7.80]`.
- Lines 30: computes `parameters.d_ambig` using `parameters.d_ambig={[1 2 3 4],[5 6],[7 8 9]}`.
- Lines 33: computes `parameters.p_shifts` using `parameters.p_shifts=[6.5 -56.9 -21.7 -23.0 11.4 54.6 18.6 16.4 17.8]`.
- Lines 36: computes `parameters.p_ambig` using `parameters.p_ambig={[3 4],[7 8 9]}`.
- Lines 39: computes `[~,~,pcs_theo,pcs_expt,~,~]` using `[~,~,pcs_theo,pcs_expt,~,~]=pcs_combi_fit(parameters)`.

## Implementation structure

- Extracting the susceptibility tensor from DFT hyperfine tensors and
- experimental paramagnetic shifts. A combinatorial procedure is used
- that cycles through ambiguous assignments.
- Calculation time: minutes.
- Read DFT HFCs in Gauss
- Isotope list
- Spin groups with identical PCS
- Diamagnetic shifts
- Diamagnetic shift ambiguities
- Paramagnetic shifts
- Paramagnetic shift ambiguities
- Run the combinatorial fitting

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `gparse()`, `pcs_combi_fit()`, `kfigure()`.
