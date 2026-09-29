# examples/dnp_liq/jdnp/fig_1_exch_and_field_scan.m

- MATLAB implementation: [examples/dnp_liq/jdnp/fig_1_exch_and_field_scan.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/dnp_liq/jdnp/fig_1_exch_and_field_scan.m)

## What it computes

This is a finite-time JDNP proton-polarisation map over external field and inter-electron exchange coupling, associated by the source with Fig. 1 in DOI: [10.1039/d1cp04186j](https://doi.org/10.1039/d1cp04186j). The source labels the calculation as hours with line-by-line plotting; this is an estimate from the source, not a measured runtime here.

## Required setup and model defaults

Run in MATLAB with Spinach and the sibling `examples/dnp_liq/jdnp/system_specification.m` on the path. The no-argument function loads `[sys,inter,bas,parameters]` from that helper. Its model is `{'1H','E','E'}`: proton chemical-shift matrix `diag([5 10 20])`, both electron g matrices `diag([2.0032 2.0032 2.0026])`, coordinates `[-3.00 0.50 1.30]`, `[0 0 -9.37]`, `[0 0 9.37]`, and initially empty scalar coupling. The helper sets Redfield plus SRFK relaxation, `equilibrium='dibari'`, `rlx_keep='labframe'`, temperature value `298`, `tau_c={500e-12}`, `srfk_tau_c={[1.0 1e-12]}`, and `srfk_mdepth{2,3}=3e9`; the units of the coordinate entries and these helper parameters are not explicitly annotated there. It uses `sphten-liouv` with no basis approximation and relaxation-integration tolerance `1e-10`; hygiene checks are disabled and output is set to hush. Reference g values are `2.00231930436256` and the mean electron g value.

The scan overrides the field at each outer-loop step and sets the {2,3} scalar coupling for each inner-loop point. Microwave power is `(2*pi*1e6)/2` rad/s; pulse duration is `20e-3` seconds.

## Grid, dynamics, and output

The field grid is 256 points from 0.25 to 3.0 T; exchange grid is 256 points from `-100e9` to `+100e9` Hz. For each field the code computes free-electron and trityl frequencies and sets the microwave offset to their difference. The exchange points run in `parfor`: build the system/basis, obtain electron `Lx/2` and `Lz` operators, equilibrium state, and proton `Lz` detection state; form the ESR Hamiltonian and relaxation superoperator; add the microwave and offset terms; then evolve `rho_eq` under `H+1i*R` for 20 ms and take the final state. The plotted value is `real(Nz'*rho)/real(Nz'*rho_eq)`. The heat map uses exchange in GHz horizontally and field in Tesla vertically, and updates after each field row. It evaluates 256×256 points in total and does not save the result matrix or figure to a file.

**Source-specific clarification:** despite the figure's “matching condition” framing, this script evaluates a fixed 20 ms endpoint polarisation normalised to equilibrium—not a steady-state solution; its progress plot is refreshed once per field row.

**Caveats:** use a MATLAB environment that supports `parfor`. Source estimates hours of computation. The map, output values, and article figure agreement were the helper's unlabelled parameter units are intentionally not guessed.