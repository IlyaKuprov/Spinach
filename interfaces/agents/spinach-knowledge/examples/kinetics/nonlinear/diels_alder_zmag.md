# examples/kinetics/nonlinear/diels_alder_zmag.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/kinetics/nonlinear/diels_alder_zmag.m`
- Signature: `diels_alder_zmag()`
- Total lines: 170

## Purpose

Time-domain Z magnetisation dynamics in the Diels-Alder cycloaddition of acetylene to butadiene, demonstrating the non-linear kinetics module. Calculation time: minutes.

## Physical / mathematical content

- Chemical-kinetics examples. The files couple spin dynamics to exchange, pumping, or nonlinear reaction networks represented by kinetic generators in Liouville space.
- Orientation or trajectory averaging is performed numerically, so grid design, weights, and integration error control matter directly to accuracy and runtime.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- Numerical integration over angles or geometry is part of the implementation, so point placement and weights are as important as the local Hamiltonian calculations.

## Implementation structure

- Time-domain Z magnetisation dynamics in the Diels-Alder cycloaddition
- of acetylene to butadiene, demonstrating the non-linear kinetics module.
- Calculation time: minutes.
- DFT import options
- Load and display acetylene (substance A)
- Load and display butadiene (substance B)
- Load and display cyclohexadiene (substance C)
- Add natural abundance ethanol (substance D)
- Merge the spin systems
- Magnet field
- Chemical parts and unit concentrations
- Basis set

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `gparse()`, `g2spinach()`, `kfigure()`, `scale_figure()`, `subplot()`, `cst_display()`, `camorbit()`, `ktitle()`, `num2cell()`, `merge_inp()`, `create()`, `basis()`, `step()`, `griddedInterpolant()`, `kxlabel()`, `kylabel()`.
