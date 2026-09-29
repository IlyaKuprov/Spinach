# examples/dnp_liq/jdnp/fig_3_time_dep_top_row.m

- MATLAB implementation: [examples/dnp_liq/jdnp/fig_3_time_dep_top_row.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/dnp_liq/jdnp/fig_3_time_dep_top_row.m)

## What it calculates

This zero-argument example is the “top row” control paired with fig_3_time_dep_bot_row.m. Its source describes the intended demonstration as removal of the second electron to show the JDNP contribution absent from the corresponding two-electron case. It plots proton polarisation time traces at three fields and cites [DOI: 10.1039/d1cp04186j](https://doi.org/10.1039/d1cp04186j). That stated purpose is not a numerical result reproduced here.

## Model edits and run setup

The function begins with system_specification() and assumes its original three-site ordering: proton at site 1, retained electron at site 2, and the electron to be removed at site 3. It keeps sites 1 and 2 in sys.isotopes, inter.zeeman.matrix, and inter.coordinates; resets inter.coupling.scalar to an empty 2-by-2 cell array; removes srfk_tau_c and srfk_mdepth; and sets inter.relaxation={'redfield'}. Thus this is not the bottom-row system with an electron merely hidden from the plot: the model data are explicitly reduced before create().

Run with MATLAB and Spinach, the JDNP helper, and Spinach plotting routines available on the path. It sets parameters.mw_pwr=2*pi*250e3, t_step=1e-3, and nsteps=300; the field grid is [0.034 0.34 3.4] T. For each field it derives mw_off from g_ref and g_trityl, then builds the ESR Hamiltonian and relaxation operator and evolves proton Lz from isotropic thermal equilibrium. Unlike the bottom-row source, it does not assign a field-matched {2,3} scalar coupling after removing the second electron.

## Output and limits

It opens a 1-by-3 figure with time in seconds and one panel per field; the time grid runs from zero through t_step*nsteps (301 samples). It plots normalised real proton Lz expectation values and does not save/export data. The helper supplies the starting system and parameters; the source does not list all helper defaults. The file's intended “vanishing JDNP” comparison has not been independently validated here, and mw_pwr has no unit stated in the source.
