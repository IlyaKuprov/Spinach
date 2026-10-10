# examples/dnp_liq/jdnp/fig_4_state_amplitudes.m

- MATLAB implementation: [examples/dnp_liq/jdnp/fig_4_state_amplitudes.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/dnp_liq/jdnp/fig_4_state_amplitudes.m)

## What it calculates

This example follows selected operator-state amplitudes for the proton plus two-electron JDNP system, using the electron singlet/triplet combinations and the nuclear-spin-resolved components defined explicitly in the source. The stated purpose is to illustrate imbalance between singlet–nuclear-alpha and singlet–nuclear-beta subspaces and the associated transient nuclear polarisation. It cites [DOI: 10.1039/d1cp04186j](https://doi.org/10.1039/d1cp04186j).

## Running assumptions and settings

The zero-argument function reads sys, inter, bas, and parameters from system_specification() and assumes proton at site 1 and electrons at sites 2 and 3. Run with MATLAB, Spinach, that JDNP helper, and Spinach plotting routines available on the path. It sets the field to 14.08 T, parameters.mw_pwr=2*pi*250e3, t_step=1e-3, and nsteps=700; the microwave offset is calculated from g_ref and g_trityl. As in the field-sweep example, the scalar coupling {2,3} is set to the sum of the isotropic electron and proton Zeeman terms. Propagation uses the ESR Hamiltonian plus relaxation and requests the multichannel evolution mode.

The 15 detection channels are ordered as [Tpa,Tpb,Tma,Tmb,T0a,T0b,Sa,Sb,SNz,TpNz,T0Nz,TmNz,E1z,E2z,Nz]. Here the first six are the electron triplet components resolved by nuclear alpha/beta, channels 7–8 are the alpha/beta singlets, 9–12 are nuclear-z-weighted singlet/triplet components, and 13–15 are electron 1, electron 2, and nuclear Lz components, respectively.

## Output and limits

The function creates a three-panel figure: panel 1 plots triplet channels 1, 3, 5 (alpha) and 2, 4, 6 (beta); panel 2 plots singlet channels 7–8; panel 3 plots channel 15 (Nz). Each trace uses the real part of the corresponding returned channel. All panels plot t_axis(2:end) against answer(:,2:end), so the initial time sample is deliberately not displayed. The time grid spans zero through t_step*nsteps (701 samples), with axes labelled in seconds. No output file is written.
