# examples/nmr_solids/case_studies/mas_powder_dd_nqi.m

Source: [mas_powder_dd_nqi.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_solids/case_studies/mas_powder_dd_nqi.m)

## Spin system

The header describes a powder MAS spectrum of a dipole-coupled pair of quadrupolar nuclei and credits the parameters to Jeongjae Lee. The system contains ²³Na and ¹⁷O at a field of 9.4 T, with coordinates (2.602, 8.750, 3.651) and (4.401, 10.184, 4.371). It supplies separate CAStep-to-NQI tensors: for ²³Na, the matrix rows are (−0.0497, 0.0520, −0.0019), (0.0520, 0.0315, 0.0027), and (−0.0019, 0.0027, 0.0182), followed by +0.1040 and spin 3/2; for ¹⁷O, the rows are (0.1580, 0.0340, −0.5562), (0.0340, −0.6005, 0.0586), and (−0.5562, 0.0586, 0.4425), followed by −0.0258 and spin 5/2. The source does not state units for the magnet setting, coordinates, or tensor arguments.

The basis is sphten-liouv with no approximation and projection +1. A possible GPU setting is present only as a commented-out line, so this script does not enable it.

## MAS acquisition and spectrum

The acquisition parameters include rate 100000, axis [√(2/3), 0, √(1/3)], maximum rank 30, a 200-point spherical grid, sweep 5 × 10⁶, 1,024 points, zero-fill to 4,096, and offset zero. The source does not state a unit for rate, sweep, or offset; it explicitly sets the plotted axis units to MHz. It selects the ¹⁷O L+ state for both initial state and receiver coil, with no decoupling list and no RF or CP pulse sequence specified.

Spinach performs single-rotor acquisition, applies exponential apodisation with parameter 6, Fourier transforms the signal, and plots the real spectrum. This is a simulation, not an experimental spectrum supplied by the source. The pair is dipole-coupled through its coordinates: `create()` invokes `dipolar()` to generate the ²³Na–¹⁷O point-dipolar tensor automatically. The explicit diagonal NQI assignments are separate; adding another off-diagonal dipolar tensor would double-count the interaction.

Zero track elimination is explicitly enabled with `zte` in `sys.enable`.
