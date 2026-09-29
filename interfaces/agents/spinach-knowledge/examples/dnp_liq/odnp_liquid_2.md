# examples/dnp_liq/odnp_liquid_2.m

- MATLAB implementation: [examples/dnp_liq/odnp_liquid_2.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/dnp_liq/odnp_liquid_2.m)

- Signature: `odnp_liquid_2()`
- Calculation time: seconds
- [MATLAB source](../../../../../examples/dnp_liq/odnp_liquid_2.m)

## Purpose

Models liquid-phase Overhauser DNP at room temperature after a perfect inversion of the electron ESR signal. Redfield relaxation represents dipolar electron–nuclear cross-relaxation. The simulation follows longitudinal signals after the inversion rather than applying a continuous microwave drive.

## Spin system and relaxation

The three spins are `1H`, `1H`, and `E`, with `sys.magnet=3.4`. The Zeeman matrices are isotropic for both protons (diagonal entries 5) and specify electron entries 2.0023, 2.0025, and 2.0027. Coordinates, in the source's stated Angstrom units, are (0, 0, 0), (0, 2, 0), and (0, 0, 1.5), respectively. The basis is `sphten-liouv` with no approximation. Relaxation uses `redfield`, Di Bari equilibrium, secular retention, temperature 298, and a 10 ps correlation time.

## Preparation and time evolution

The script creates the spin system and basis, obtains `rho_eq=equilibrium(spin_system)`, then applies `step(spin_system,Lx,rho_eq,pi)` using the electron `Lx` operator. This prepares `rho0` by a pi-radian inversion. The ESR-context call is `liquid(spin_system,@dnp_time_dep,parameters,'esr')`; the parameter structure selects the electron, supplies `rho0`, and defines proton 1, proton 2, and electron `Lz` states as the three coil channels. Microwave power and offset are both zero. The time step is `1e-6` and the number of steps is `1e3`.

## Output and scope

The first plot is the real electron longitudinal signal (the third coil channel); the second plots the two proton longitudinal signals. The plotted abscissa is `linspace(0,1000,1001)`, labelled in microseconds, and the ordinates are longitudinal `Lz` signals. This is a fixed three-spin geometry with one correlation time and no microwave drive; it is not a field or frequency sweep.
