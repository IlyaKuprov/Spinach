# examples/dnp_liq/jdnp/fig_2_tau_and_field_traject.m

- MATLAB implementation: [examples/dnp_liq/jdnp/fig_2_tau_and_field_traject.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/dnp_liq/jdnp/fig_2_tau_and_field_traject.m)

## What it calculates

This zero-argument example constructs a liquid-state JDNP time trace for six static-field settings, with four rotational correlation times at each field. The inter-electron scalar coupling is reset at each field to the sum of the isotropic electron Zeeman frequency and the proton Zeeman frequency; the resulting figure compares time-dependent proton polarisation across fields and correlation times. The source points to [DOI: 10.1039/d1cp04186j](https://doi.org/10.1039/d1cp04186j).

## Running assumptions and settings

Run with MATLAB and Spinach available on the path, and with the JDNP example helper system_specification() and Spinach plotting helpers accessible. The function takes no arguments; it obtains sys, inter, bas, and parameters from that helper. Its indexing assumes the proton is site 1 and the two electrons are sites 2 and 3, including inter.zeeman.matrix{2} and scalar coupling {2,3}.

The source sets parameters.mw_pwr=2*pi*250e3, t_step=1e-3, and nsteps=200. The six fields are [0.5 3.4 7.0 11.7 14.1 23.5] T. tau_c is [300 400 500 600] ps (stored in seconds as [300e-12 400e-12 500e-12 600e-12]). At each field, the microwave offset is formed from g_ref and g_trityl; the code then loops over tau_c, builds the Spinach system/basis, and propagates with the ESR Hamiltonian plus relaxation and microwave terms. The proton Lz expectation is made real and divided by its thermal-equilibrium value.

## Output and limits

The routine opens a 2-by-3 figure, one panel per field, with four correlation-time traces in each. It constructs 201 times from zero through t_step*nsteps, displays them as milliseconds (1e3*t_axis), and fixes the plotted vertical range to [-250, 50]. It creates a figure only; it does not save the trace or export the plot. The numeric setting mw_pwr is recorded as written because the source does not state its unit. No simulation or numerical result is claimed to have been reproduced here.
