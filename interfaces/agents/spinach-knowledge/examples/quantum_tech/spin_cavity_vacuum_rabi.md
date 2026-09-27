# examples/quantum_tech/spin_cavity_vacuum_rabi.m

- Signature: `spin_cavity_vacuum_rabi()`

## Purpose

Vacuum Rabi oscillation between an electron spin and a microwave cavity mode in the Jaynes–Cummings approximation. This is the one-spin limit of the spin-ensemble cavity experiments of Schuster et al. and Kubo et al., Phys. Rev. Lett. 105, 140501 and 140502 (2010). Calculation time: seconds.

## Physical / mathematical content

- A resonant cavity exchanges a single excitation with an electron spin. The example tracks spin and cavity excitation populations and checks both visible transfer and conservation of population in the active doublet.

## Numerical / algorithmic content

- Spinach builds a `zeeman-hilb` model with no basis approximation and propagates the initial spin excitation through the cavity device context. The sequence uses 501 points, corresponding to the plotted 0–500 ns interval.

## Implementation structure

- The source uses isotopes `{'E','C3'}`, a resonant cavity mode at zero rotating-frame frequency, and exchange coupling `8e6`. The initial state is `{'ZL2','BL1'}` on the spin and cavity. Projectors `{'ZL2','E'}` and `{'ZL1','BL2'}` measure the spin and cavity populations, respectively.
