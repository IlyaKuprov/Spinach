# examples/dnp_liq/ccdnp/freq_scan_main_text.m

- MATLAB implementation: [examples/dnp_liq/ccdnp/freq_scan_main_text.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/dnp_liq/ccdnp/freq_scan_main_text.m)

freq_scan_main_text() computes a steady-state proton DNP response over microwave-frequency offset and magnetic field for a liquid three-spin system (1H, E, E). The source describes two exchange-coupled electrons, each dipolar-coupled to the nucleus, and cites [J. Magn. Reson. 2021, 106940](https://doi.org/10.1016/j.jmr.2021.106940). It requires Spinach, MATLAB, the dnp_freq_scan callback on the path, and Parallel Computing Toolbox for its parfor; there are no external data files.

The model uses proton Zeeman eigenvalues [0 10 20], electron g-tensor eigenvalues [2.0034 2.0038 2.0038] and [2.0057 2.0030 2.0030] with the source's stated Euler triples [-0.872 -0.013 0.868] and [-1.145 0.061 1.143]. Scalar electron exchange is set to 3e6; the three coordinate rows are [0 0 0], [5.090 0.010 0.958], and [-5.090 0.061 1.032]. These numerical interaction and coordinate values are given without units in the source. The basis is sphten-liouv with no approximation. Relaxation is Redfield with zero equilibrium, secular retention, temperature 298, correlation time 100e-12, and integration tolerance 1e-10 (the source comments that this tolerance needs to be this tight).

Sequence settings are electron irradiation, mw_pwr=2*pi*1e6, method lvn-backs, and needs={'rho_eq'}; g_ref is the mean of the first electron's g eigenvalues. The frequency-offset grid is 2*pi*linspace(-5,10,512)*1e6; the field grid is linspace(1,20,128) Tesla. A parfor loop builds the system and steady state at each field, calls liquid(...,@dnp_freq_scan,...,'esr'), then divides each frequency trace by the equilibrium coil projection coil'*rho_eq.

The output is a figure of real(answer) against field and frequency offset, with a colorbar labelled steady-state 1H DNP; the plotted offset is converted to MHz relative to the isotropic g reference. No result array is saved. The source comments “Calculation time: seconds”; this is not a reproduced runtime or a performance guarantee.
