# examples/dnp_liq/jdnp/fig_5_spatial_distribution.m

- MATLAB implementation: [examples/dnp_liq/jdnp/fig_5_spatial_distribution.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/dnp_liq/jdnp/fig_5_spatial_distribution.m)

- Signature: fig_5_spatial_distribution()
- Reference: [Concilio et al., *Physical Chemistry Chemical Physics* (2022)](https://doi.org/10.1039/d1cp04186j)
- The source comment gives a calculation time of seconds.

## Purpose and run context

Maps the proton longitudinal-polarisation signal around the radical pair for the JDNP parameter choice described by the source. Call fig_5_spatial_distribution(); it obtains sys, inter, bas, and parameters from system_specification(), then creates two position maps. The source introduction frames this as illustrating JDNP persistence under liquid-phase position and orientation averaging; the computation here explicitly scans proton position in two planes.

## Model and settings

The routine sets sys.magnet=14.08, microwave power to 2*pi*250e3, and evolution time to 20e-3 (identified as 20 ms in the source header). The microwave offset is 2*pi*(f_trityl-f_free), with both frequencies obtained using g2freq and the helper's reference g-factors. It sets the electron-pair scalar interaction at inter.coupling.scalar{2,3} (the coupling between electron spins 2 and 3) to the sum of the electron isotropic Zeeman frequency and proton Zeeman frequency. The shared helper supplies a three-spin 1H,E,E system and a complete sphten-liouv basis.

Each coordinate array is linspace(-30,30,30). The source labels coordinate axes in Angstrom. The first scan places the proton at [X(n),Y(k),0] (the Z=0 plane); the second uses [X(n),0,Z(k)] (the Y=0 plane). Thus the base proton coordinate from system_specification() is replaced at each grid point. Both inner scans use parfor.

## Propagation and output

At each point the code creates the spin system and basis, forms electron Lx and Lz operators and proton Lz detection, obtains isotropic thermal equilibrium, and builds the ESR Hamiltonian and relaxation superoperator. It adds mw_pwr*Lx + mw_off*Lz for the electron and propagates rho_eq with H+1i*R for the selected evolution time. The map value is real(Nz'*rho)/real(Nz'*rho_eq).

The figure contains colour maps for the Z=0 and Y=0 planes; the colour-bar label identifies proton DNP at 20 ms. The XZ-panel colour limits are set to [-200,1]; no explicit colour limits are assigned to the XY panel. The source scans position; it contains no separate explicit orientation-sampling loop.
