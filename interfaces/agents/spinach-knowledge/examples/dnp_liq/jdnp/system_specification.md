# examples/dnp_liq/jdnp/system_specification.m

- MATLAB implementation: [examples/dnp_liq/jdnp/system_specification.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/dnp_liq/jdnp/system_specification.m)

- Signature: [sys,inter,bas,parameters]=system_specification()
- Reference: [Concilio et al., *Physical Chemistry Chemical Physics* (2022)](https://doi.org/10.1039/d1cp04186j)

## Purpose and use

Returns the shared three-spin model used by the adjacent JDNP examples: one proton and two electrons, their interactions and relaxation settings, a complete basis specification, and reference g-factors. The scalar-coupling container is initialised but left for the calling simulation to set. fig_5_spatial_distribution() requests all four outputs; fig_6_microwave_free() requests the first three and then sets its own scalar coupling and correlation time.

## Spin system and interactions

The isotopes are {'1H','E','E'}. The proton Zeeman tensor is diag([5 10 20]); each electron g-tensor is diag([2.0032 2.0032 2.0026]). Coordinates are assigned as proton [-3.00 0.50 1.30], electron 2 [0.00 0.00 -9.37], and electron 3 [0.00 0.00 +9.37]. The scalar-coupling field is a 3-by-3 cell array without assigned values.

## Relaxation, basis, and returned parameters

The relaxation selection is {'redfield','SRFK'}, the equilibrium method is dibari, and relaxation is retained in the lab frame. The source sets temperature to 298, tau_c={500e-12}, srfk_tau_c={[1.0 1e-12]}, and srfk_mdepth{2,3}=3e9 in a 3-by-3 cell array. The basis uses sphten-liouv formalism with approximation='none'; sys.tols.rlx_integration is 1e-10. It disables hygiene and sets output to hush.

The returned parameters.g_ref is 2.00231930436256, and parameters.g_trityl is the mean of the diagonal entries of electron 2's g-tensor. No units for the tensor entries, temperature value, or correlation-time settings are stated in this source.
