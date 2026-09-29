# examples/dnp_liq/jdnp/fig_6_microwave_free.m

- MATLAB implementation: [examples/dnp_liq/jdnp/fig_6_microwave_free.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/dnp_liq/jdnp/fig_6_microwave_free.m)

- Signature: fig_6_microwave_free()
- The source comment gives a calculation time of minutes.
- The shared model is described in system_specification.m, which cites [Concilio et al., *Physical Chemistry Chemical Physics* (2022)](https://doi.org/10.1039/d1cp04186j).

## Purpose and run context

Demonstrates a microwave-free JDNP field-ramp trajectory. Call fig_6_microwave_free(); it obtains the spin system, interactions, and basis from system_specification(), modifies the electron-pair scalar interaction and correlation-time setting, and propagates the initial thermal-equilibrium state through the ramp.

## Field ramp and propagation

The start, match, and final field values are 14.09, 11.74, and 9.39. The scalar interaction between electron spins 2 and 3 is set to match_field*(spin('E')+spin('1H'))/(2*pi); the correlation-time entry is set to 2.2e-9. The system is constructed at the start field, its basis is built, and the lab-frame Hamiltonian is used to obtain the initial equilibrium state.

The field sequence is linspace(start_field,final_field,211), with dt=1e-4. The initial state plus 211 propagated states form the trajectory. At each field value the code updates sys.magnet, recreates the system and basis, constructs the lab-frame Hamiltonian and relaxation superoperator, then calls evolution for one step using H+1i*R. No microwave drive is added. The plotted time axis is labelled in seconds and spans 0 to 0.0211.

## Observables and limitation

The code constructs product-state operators for the alpha- and beta-manifold singlet/triplet components, additional nuclear-N_z-resolved singlet/triplet components, and the longitudinal operators of both electrons and the nucleus. It computes their real projections onto the trajectory. Three panels display the alpha/beta triplet populations, the alpha/beta singlet populations, and nuclear N_z, respectively.

The source header states that the field ramp and unequal singlet-alpha/singlet-beta relaxation rates produce enhancement beyond the Boltzmann level at both fields. The plotted nuclear panel is the N_z trajectory; the script does not plot an explicit Boltzmann reference curve.
