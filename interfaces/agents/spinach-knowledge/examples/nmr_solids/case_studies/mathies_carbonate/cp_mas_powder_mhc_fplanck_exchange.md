# examples/nmr_solids/case_studies/mathies_carbonate/cp_mas_powder_mhc_fplanck_exchange.m

- Signature: `cp_mas_powder_mhc_fplanck_exchange()`
- Source: https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_solids/case_studies/mathies_carbonate/cp_mas_powder_mhc_fplanck_exchange.m

## Purpose

Simulates 13C-detected cross-polarisation (CP) contact curves for H1, H4, and C19 sites in monohydrocalcite under magic-angle spinning (MAS), with chemical exchange between two three-spin endpoint configurations. The source cites https://doi.org/10.1038/s41467-023-44381-x and estimates hours of CPU time, with a GPU much faster.

## Spin and exchange model

The calculation reads the CASTEP-derived `mhc.magres` file and drops O and Ca. It constructs six spins: endpoint 1 has H1, H4, and C19 (indices 1-3), and endpoint 2 has H4, H1, and C19 (indices 4-6). The isotope sequence is 1H, 1H, 13C in each endpoint; the repeated Cartesian coordinates and shielding tensors select the corresponding sites from the input. Proton tensors use the source's Huang et al. ACIE 2021 parametrisation, `29.25*eye(3)-cst`, and the C19 tensors use `169.86*eye(3)-cst`. The two kinetic parts are `[1 2 3]` and `[4 5 6]`, with initial concentrations `[1 1]` in arbitrary units.

It runs six exchange-rate cases: 10, 100, 1,000, 10,000, 100,000, and 1,000,000 Hz. The source labels the spectrometer setting as 400 MHz NMR and sets `sys.magnet=9.4`.

## MAS, CP contact, and detection

The source sets the MAS rate parameter to 10,000 (no unit is written beside this assignment) about axis `[1 1 1]`; the powder grid is `rep_2ang_800pts_sph`, with maximum rank 7. The offset values are `[2,000 10,000]`; the source comments these as 5 ppm for 1H and 100 ppm for 13C. The listed high-power field is 83,000 Hz and the two CP powers are `[60,000 50,000]` Hz. The source does not assign those two CP entries to individual channels in a comment.

For each exchange rate, the code calls `singlerot` with `@cp_contact_soft`, using 1,000 steps of 10 microseconds (a 10 ms contact-time span); it also requests `parameters.needs={'iso_eq'}`. This is the CP contact simulation selected by the source, rather than a separate hand-coded pulse train in the example. It detects with `coil_state(spin_system,'L+','13C','exact')`; this is the 13C transverse single-quantum observable. The plotted output is the real simulated 13C signal in arbitrary units against contact time, with one curve per exchange rate.

## Inputs and outputs

The input is the CASTEP-derived structural/shielding data in `mhc.magres`; this source does not load an experimental FID or spectrum. The output is a simulated contact-curve plot. The GPU enable line is commented out, and the source notes that a GPU can accelerate the hours-long calculation.

Zero track elimination is explicitly enabled with `zte` in `sys.enable`.

The rate matrix is represented by two explicit first-order reaction records. Atom matching pairs equal-position spin indices in the two declared parts in both directions; the permuted geometry/tensors represent the conformational exchange. Detection uses unweighted `coil_state`.
