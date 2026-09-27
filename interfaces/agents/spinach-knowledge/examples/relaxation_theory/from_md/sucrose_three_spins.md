# examples/relaxation_theory/from_md/sucrose_three_spins.m

- Signature: `sucrose_three_spins()`

## Purpose

One of the calculations reported in the JMR paper with Jim Prestegard: a three-spin subsystem from the glucose ring of sucrose, illustrating the incorrect viscosity of TIP3P water. The experimental value of tau_c for sucrose is around 90 ps (correctly reproduced by OPC and TIP5P water), but TIP3P only agrees with Redfield theory when tau_c is set to 37 ps in the latter. Here, the numerical relaxation superoperator computed from a long MD trajectory is compared with the analytical one computed using the isotropic rotational diffusion approximation.
## Implementation structure

- One of the calculations reported in the JMR paper with Jim Prestegard:
- a three-spin subsystem from the glucose ring of sucrose -an illustra-
- tion of incorrect viscosity of TIP3P water.
- The experimental value of tau_c for sucrose is around 90 ps (and this
- is correctly reproduced by OPC and TIP5P water), but TIP3P only agrees
- with Redfield theory when tau_c is set to 37 ps in the latter.
- Here, numerical relaxation superoperator computed from a long MD tra-
- jectory is compared with the analytical one computed using the isotro-
- pic rotational diffusion approximation.
- Calculation time: minutes, with most of the time spent
- computing MD frame Hamiltonians
- Three-spin system
