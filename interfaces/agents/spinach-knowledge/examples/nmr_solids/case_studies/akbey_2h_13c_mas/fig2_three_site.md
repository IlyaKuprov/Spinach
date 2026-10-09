# examples/nmr_solids/case_studies/akbey_2h_13c_mas/fig2_three_site.m

- MATLAB implementation: [examples/nmr_solids/case_studies/akbey_2h_13c_mas/fig2_three_site.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_solids/case_studies/akbey_2h_13c_mas/fig2_three_site.m)

Source: https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_solids/case_studies/akbey_2h_13c_mas/fig2_three_site.m

This example simulates the three-site position exchange of a deuterium nucleus. Its comments describe the sites as differing in chemical shift and quadrupolar-tensor orientation and state that the calculation is intended to reproduce Figure 2 of the cited paper. The spin-system entries are three ²H chemical sites, with equal concentrations and an exchange-rate matrix formed by multiplying the stated matrix by 1e4. The driver constructs a separate quadrupolar tensor for each site using weblab2nqi(0.16e6, 0, 1, 0, acos(1/3), 2*pi/3). It sets the chemical exchange-rate matrix to 1e4 times [-2, 1, 1; 1, -2, 1; 1, 1, -2] and concentrations to [1, 1, 1]. The source does not annotate units for these inputs. It sets the magnetic-field parameter to 9.4; the source does not annotate a unit for this value, the quadrupolar-helper inputs, or the exchange rates.

The experiment is a simulated one-dimensional rotor-synchronised ²H acquisition, not processing of experimental data. It starts and detects ²H transverse coherence using L+ operators, then calls the single-rotor acquisition routine. The rotor parameter is 8500 and the rotor axis is [1 1 1]; their units are not stated in the source. The powder grid is rep_2ang_800pts_sph and the maximum rank is 25. The spectral setup uses offset 0, a 0.4e6 sweep input, 1024 acquired and zero-filled points, and an axis labelled in kHz; the axis is inverted. The source does not state units for the sweep or rotor-rate inputs. No explicit RF pulse is specified in this driver.

No apodisation is applied. The output plot is the real part of the Fourier-transformed simulated FID. The source cites Figure 2 of Umit Akbey et al., Journal of Magnetic Resonance (2021), DOI: https://doi.org/10.1016/j.jmr.2021.106974.
