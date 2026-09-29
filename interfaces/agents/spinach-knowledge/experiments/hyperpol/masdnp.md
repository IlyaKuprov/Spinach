# experiments/hyperpol/masdnp.m

- MATLAB implementation: [experiments/hyperpol/masdnp.m](https://github.com/IlyaKuprov/Spinach/blob/main/experiments/hyperpol/masdnp.m)

- Signature: `dnp=masdnp(spin_system,parameters)`

## Purpose and physical scope

This routine simulates magic-angle-spinning DNP and returns a powder-orientation-weighted enhancement ratio. For each orientation it constructs an ESR rotor stack, applies microwave excitation and relaxation through a rotor period, forms an effective period generator, evolves the thermal-equilibrium state for the specified microwave duration, then averages the detected magnetisation over a rotor period and divides by the corresponding equilibrium `coil` signal. The spin-system model supplies the couplings; the routine does not construct hyperfine tensors or run ESEEM/ENDOR or image reconstruction.

## Inputs

- `parameters.spins`: spin labels for microwave irradiation; the implementation uses the first label, `parameters.spins{1}`, for the microwave frequency and operator.
- `parameters.rate`: MAS rate in Hz; `parameters.axis`: spinning-axis direction vector; `parameters.max_rank`: integer rotor-discretisation rank.
- `parameters.mw_pwr`: microwave power in radians per second; `parameters.mw_frq`: microwave frequency in Hz; `parameters.mw_time`: irradiation/equilibration duration in seconds.
- `parameters.grid`: name of a spherical averaging grid available in the Spinach grids folder; `parameters.coil`: detection state; `parameters.verbose`: diagnostic-output flag (0 or 1).
- `spin_system` must use `sphten-liouv` or `zeeman-liouv` formalism. Call `masdnp` directly, not through a context wrapper.

## Calculation and output

The routine builds the microwave operator from the first selected spin's `L+` and `L-` operators, constructs relaxation and a lab-frame thermal-equilibrium state, and iterates over the weighted spherical grid. The rotor propagator is assembled from step propagators over one rotor period; its logarithm supplies an effective generator for `evolution` over `parameters.mw_time`. The function then averages the detection-state projection over one rotor period and accumulates the orientation-weighted ratio `Hz_dnp/Hz_eq`. The output dnp is a scalar enhancement ratio, not an absolute polarisation or a time/frequency spectrum.

## Model limits

The source recommends increasing rotor rank and spherical-grid size until the answer stops changing, noting that both may need to be very large. No convergence is implied by the example values below. The routine requires the named spherical-grid file to be installed under the Spinach root directory.

## Source-coded numerical example and DOI

examples/dnp_mas/solid_effect_mas_powder.m sets parameters.rate=12.5e3 Hz, parameters.max_rank=800, parameters.mw_pwr=2*pi*0.85e6 radians per second, parameters.mw_frq=-263.366e9 Hz and parameters.mw_time=1.0 seconds, using grid 'rep_2ang_100pts_sph'. The example comments identify the simulation as based on Fred Mentink-Vigier's paper and explicitly note that Spinach rotation conventions differ. The cited DOI is <https://doi.org/10.1016/j.jmr.2015.07.001>. These are example inputs, not experimental measurements or a claimed reproduced result.

## Source and attribution

- Source: `experiments/hyperpol/masdnp.m`
- <https://spindynamics.org/wiki/index.php?title=masdnp.m>
- Source attributions: ilya.kuprov@weizmann.ac.il; fmentink@magnet.fsu.edu
