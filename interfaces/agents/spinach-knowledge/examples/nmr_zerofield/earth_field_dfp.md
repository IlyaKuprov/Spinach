# examples/nmr_zerofield/earth_field_dfp.m

- Signature: `earth_field_dfp()`
- Source: [`examples/nmr_zerofield/earth_field_dfp.m`](https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_zerofield/earth_field_dfp.m)

## Experiment and spin model

This example simulates Earth-field NMR for 2,6-difluoropyridine and reproduces the simulated spectra in Figure 7 of [the cited paper](https://doi.org/10.1016/j.jmr.2023.107540), omitting the weighted addition of the uncoupled 1H signal. It is not a zero-field experiment or an imported measured spectrum. The source credits Adam Altenhof and Derrick Kaseman. The field is set so the proton frequency is 2312.45 Hz. The six sites are H3, H4, H5, F2, F6 and N1, with isotopes 1H, 1H, 1H, 19F, 19F and 14N. Chemical-shift entries are 6.98, 8.06, 6.98, −70.69, −70.69 and 0 ppm, respectively.

The source labels its scalar couplings experimental and encodes (in Hz): F2–N1 and F6–N1, 37.35 each; H3–N1 and H5–N1, each set to the mean (1.03 + 0.71)/2 = 0.87; H3–F2 and H5–F6, −2.47 each; H4–F2 and H4–F6, 8.08 each; H5–F2 and H3–F6, 1.29 each; F2–F6, −12.23; H3–H4 and H5–H4, 7.92 each; and H3–H5, 0.55. For the two proton–nitrogen values the source comments that Table 2 “makes no sense”; the mean is the value actually implemented.

## Relaxation sweep and spectrum

The sphten-liouv basis is exact (no approximation). Relaxation uses T1/T2 rates, diagonal retention and zero equilibrium. The sequence uses a 90-degree flip, uniaxial detection, 1H channel, zero offset, 10 kHz sweep, 60,000 acquired points and zero filling to 2¹⁹ points. For each 14N relaxation rate in [0, 10, 100, 500, 1,000, 10,000, 100,000, 1,000,000] Hz, the five non-nitrogen sites receive R1 = R2 = 0.22 Hz and 14N receives the swept rate. The script calls the liquid simulator with the zero-field sequence function in the lab frame, Fourier transforms the FID, normalises the real spectrum to its maximum, offsets each trace vertically, and plots 2303–2322 Hz. The source comments that a GPU is needed and estimates seconds, but its explicit GPU-enable line is commented out; neither hardware requirement nor runtime was validated here. It does not report a field-drop schedule.
