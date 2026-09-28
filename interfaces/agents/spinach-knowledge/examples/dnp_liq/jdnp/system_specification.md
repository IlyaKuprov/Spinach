# examples/dnp_liq/jdnp/system_specification.m

- Signature: `[sys,inter,bas,parameters]=system_specification()`
- Reference: [Concilio et al., *Phys. Chem. Chem. Phys.* (2022)](https://doi.org/10.1039/d1cp04186j)

## Purpose

Defines the shared three-spin (one proton and two electrons) model used by the adjacent JDNP examples, including its Zeeman tensors, coordinates, relaxation settings, and basis. The scalar-coupling matrix is initialized empty for the calling example to populate.

## Spin system and interactions

The isotopes are `1H, E, E`. The proton chemical-shift tensor is diagonal with entries 5, 10, and 20; both electron g-tensors are axial, with diagonal entries 2.0032, 2.0032, and 2.0026. The proton is placed at (-3.00, 0.50, 1.30) Å, while the electrons lie on the z axis at -9.37 and +9.37 Å.

## Relaxation and basis

The model uses Redfield relaxation with the SRFK mechanism, Di Bari equilibrium, lab-frame relaxation retention, temperature 298 K, and a 500 ps correlation time. The SRFK correlation-time and modulation-depth settings are specified explicitly in the source, with the latter applied to the electron pair. The basis uses the sphten-liouv formalism without an approximation; the relaxation-integration tolerance is 1e-10. Reference electron g-factors for the free electron and trityl are returned in `parameters` for use by the simulations.
