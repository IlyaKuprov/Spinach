# examples/dnp_liq/jdnp/fig_3_time_dep_bot_row.m

- MATLAB implementation: [examples/dnp_liq/jdnp/fig_3_time_dep_bot_row.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/dnp_liq/jdnp/fig_3_time_dep_bot_row.m)

## What it calculates

This is the two-electron (“bottom row”) JDNP time-trace example, intended to be read alongside fig_3_time_dep_top_row.m, which removes one electron as a control. It computes proton Lz polarisation at three fields while setting the inter-electron scalar coupling to its field-dependent matching value. The source cites [DOI: 10.1039/d1cp04186j](https://doi.org/10.1039/d1cp04186j).

## Running assumptions and settings

The zero-argument function obtains sys, inter, bas, and parameters from system_specification(). Run it with MATLAB, Spinach, the JDNP example helper, and Spinach plotting functions on the path. The source expects the helper's three-site arrangement to support proton detection on site 1 and the two electrons at sites 2 and 3 with the electron Zeeman matrix at site 2 used to set the scalar coupling {2,3}.

It sets parameters.mw_pwr=2*pi*250e3, t_step=1e-3, and nsteps=300. The field grid [0.034 0.34 3.4] is explicitly in tesla. For each field, it calculates the microwave offset from g_ref and g_trityl, and sets inter.coupling.scalar{2,3} to the isotropic electron Zeeman term plus the proton Zeeman term. It constructs the Spinach system and basis, obtains the ESR Hamiltonian and relaxation superoperator, adds the microwave terms, and evolves the proton Lz detection state from isotropic thermal equilibrium. The reported trace is the real expectation divided by coil'*rho_eq.

## Output and limits

The output is a 1-by-3 figure, with one panel per field. The time axis has 301 samples from zero through t_step*nsteps and is plotted in seconds; each panel uses ylim([-150 10]). There is no file export in the function. The scalar-coupling assignment and the 0.3-second end time follow directly from the source settings; no simulation result has been independently reproduced. mw_pwr is retained as the source expression because no unit is given there.
