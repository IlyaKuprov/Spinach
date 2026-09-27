# examples/dnp_liq/jdnp/fig_6_microwave_free.m

- Signature: `fig_6_microwave_free()`
- Calculation time: minutes

## Purpose

Simulates the microwave-free JDNP field-ramp example. The initial equilibrium state is propagated while the magnetic field is changed from 14.09 T to 9.39 T, and the script plots the singlet and triplet populations resolved by nuclear-spin projection alongside the nuclear magnetisation.

## Model and ramp

The spin system, interactions, and basis come from `system_specification()`. The electron-proton scalar coupling is set using the 11.74 T midpoint field, and the correlation time is set to 2.2 ns. The script builds the initial Spinach system at 14.09 T in the lab frame and computes its thermal-equilibrium state. It then samples a 211-point linear field grid ending at 9.39 T, with a 0.1 ms propagation step at each grid point.

At each field value, the Hamiltonian and relaxation superoperator are rebuilt and the state is advanced with `evolution`; no microwave drive is added. The trajectory is projected onto explicitly constructed singlet, triplet, electron, and nuclear-spin operators. Three panels show the alpha/beta triplet populations, the alpha/beta singlet populations, and the nuclear `N_z` signal versus the ramp time.
