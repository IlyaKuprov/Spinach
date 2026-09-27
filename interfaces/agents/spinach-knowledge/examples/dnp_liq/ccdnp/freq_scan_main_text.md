# examples/dnp_liq/ccdnp/freq_scan_main_text.m

- Signature: `freq_scan_main_text()`

## Purpose

Steady-state nuclear magnetisation as a function of microwave-frequency offset and magnetic field in a DNP experiment with two exchange-coupled electrons, both coupled to a nucleus by dipolar interactions. Further particulars: https://doi.org/10.1016/j.jmr.2021.106940. Calculation time: seconds.

## Physical / mathematical content

- Three-spin model comprising one proton and two electrons. The electron pair has 3 MHz exchange coupling; coordinates define the electron–nuclear dipolar couplings. Relaxation uses Redfield theory with a 100 ps correlation time at 298 K.

## Numerical / algorithmic content

- Computes a steady-state map over 512 microwave-frequency offsets from −5 to 10 MHz and 128 magnetic fields from 1 to 20 T. A parfor loop distributes field points, and the proton signal is normalized to the thermal-equilibrium reference.

## Implementation structure

- Steady state nuclear magnetisation as a function of microwave frequency
- offset and the magnet field in a DNP experiment with two electrons con-
- nected by exchange coupling, both coupled to a nucleus by dipolar coup-
- lings. Further particulars in:
- Calculation time: seconds
- Spin system
- Zeeman interactions
- Exchange coupling
- Coordinates
- Basis set
- Relaxation theory
- Sequence parameters
