# examples/nmr_solids/case_studies/akbey_2h_13c_mas/fig2_two_site.m

- MATLAB implementation: [examples/nmr_solids/case_studies/akbey_2h_13c_mas/fig2_two_site.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_solids/case_studies/akbey_2h_13c_mas/fig2_two_site.m)

Source: https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_solids/case_studies/akbey_2h_13c_mas/fig2_two_site.m

This example simulates two-site position exchange of a deuterium nucleus. The comments say the sites differ in chemical shift and quadrupolar-tensor orientation, and identify Figure 2 of the cited paper as the reproduction target. The model contains two ²H chemical sites with equal concentrations and two directed exchange records with population generator 1e4 times [-1, 1; 1, -1]. It constructs one quadrupolar tensor per site using weblab2nqi(0.16e6, 0, 1, 0, acos(1/3), 2*pi/3). The source does not annotate units for the helper inputs. The magnetic-field parameter is set to 9.4. The source does not state units for that value, the helper inputs, the exchange rates, or the rotor-rate input.

This is a simulated one-dimensional rotor-synchronised ²H acquisition, not a measured spectrum supplied as input. The initial and detected operators are both ²H L+ coherence; acquisition is delegated to the single-rotor routine. The rotor parameter is 8500, the rotor axis is [1 1 1], and the powder grid is rep_2ang_800pts_sph with maximum rank 25. The spectral setup uses offset 0, a 0.4e6 sweep input, 1024 acquired and zero-filled points, and an axis labelled in kHz; the axis is inverted. No explicit RF pulse is set by this driver, and units for the sweep and rotor rate are not annotated.

No apodisation is applied; the plotted signal is the real part of the Fourier-transformed simulated FID. The source cites Figure 2 of Umit Akbey et al., Journal of Magnetic Resonance (2021), DOI: https://doi.org/10.1016/j.jmr.2021.106974.

The exchange network is declared as directed first-order reaction records, one for each ordered pair of distinct sites, with rate `1e4` and a one-spin atom match. Initial coherence uses concentration-weighted `state`; detection uses unweighted `coil_state`.
