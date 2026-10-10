# examples/quantum_tech/transmon_cavity_swap.m

Source: [examples/quantum_tech/transmon_cavity_swap.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/quantum_tech/transmon_cavity_swap.m)

- Signature: `transmon_cavity_swap()`

## Model

This is a coherent Jaynes–Cummings-limit vacuum-Rabi exchange between a truncated transmon oscillator and a truncated cavity mode, not an electron-spin or defect calculation. The isotope labels `T3` and `C3` give three-level representations for the transmon and cavity. At zero magnetic field both rotating-frame mode frequencies are set to zero, so the modes are resonant; the transmon anharmonicity is `-250e6` (−250 MHz) and the exchange coupling is `20e6` (20 MHz). The source configures no external drive or dissipative terms.

The Zeeman-Hilbert basis is used without approximation. The initial state is transmon `BL2` with cavity `BL1`, i.e. one transmon excitation and the cavity in its lowest state. A cavity-context device trajectory evolves this initial condition. The trajectory settings include `sweep=2e9` and `npoints=301`; the plotted time coordinates are explicitly 0–150 ns at 301 points.

## Observable and plot

The code evaluates transmon and cavity excitation populations from separate coil operators and plots both against time. It also checks that cavity population reaches at least 0.95, transmon population falls to at most 0.05, and the two populations sum to one within `1e-6`. The curves represent ideal coherent excitation exchange in this finite model; these source-level checks do not establish measured device performance or fidelity.

The source relates the example to the circuit-QED Jaynes–Cummings model and cites Blais et al., *Reviews of Modern Physics* **93**, 025005 (2021) ([DOI](https://doi.org/10.1103/RevModPhys.93.025005)).
