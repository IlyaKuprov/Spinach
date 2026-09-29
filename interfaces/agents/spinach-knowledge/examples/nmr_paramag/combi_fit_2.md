# examples/nmr_paramag/combi_fit_2.m

- Signature: `combi_fit_2()`
- Source: [combi_fit_2.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_paramag/combi_fit_2.m)

## Purpose and scope

A combinatorial fit of DFT hyperfine-coupling (HFC) data and experimental paramagnetic shifts, with ambiguous assignments enumerated. The source describes the aim as extracting the susceptibility tensor and estimates a runtime of minutes. It does not identify a protein, residue, paramagnetic site, or metal species; this example is not source-identified as a carbonic-anhydrase case.

## Inputs and assignment data

- Reads HFC matrices from `l2_parker_funk_yb.log` using `gparse`; the source explicitly says HFCs are read in Gauss.
- Declares 27 nuclei, all isotope `1H`, grouped into nine PCS-equal sets: `[28 22 30]`, `[25 33 81]`, `[80 32 24]`, `[27 82 21]`, `[36 58 85]`, `[35 57 84]`, `[2 88 39]`, `[4 90 41]`, and `[61 92 43]`.
- Diamagnetic-shift values are `[3.62 2.65 2.65 2.86 4.95 4.10 7.40 8.00 7.80]`; assignment groups are `[1 2 3 4]`, `[5 6]`, and `[7 8 9]`.
- Paramagnetic-shift values are `[3.8 -13.7 -0.6 -3.4 5.2 20.7 10.7 10.3 10.7]`; the ambiguity group is `[3 4]`. The source labels the plotted PCS axes in ppm, but does not separately state units for these input arrays.

## Fit, output, and limits

Passes the assembled structure to `pcs_combi_fit(parameters)`. It retains the third and fourth returned values as theoretical and experimental PCS and plots them against each other, with a diagonal reference from -20 to 20. The first two and final two returned values are discarded, so this wrapper does not expose or save a fitted tensor or coordinates. No magnetic-field or temperature value, electron-coordinate input, or spectrum simulation appears in the source.
