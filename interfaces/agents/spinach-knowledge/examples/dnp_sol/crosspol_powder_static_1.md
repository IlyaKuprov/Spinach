# examples/dnp_sol/crosspol_powder_static_1.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/dnp_sol/crosspol_powder_static_1.m`
- Signature: `crosspol_powder_static_1()`
- Total lines: 53

## Purpose

E-15N cross-polarization experiment in the doubly rotating frame. Static powder simulation. Calculation time: seconds

## Physical / mathematical content

- Solid-state DNP examples. These files model microwave-driven electron-nuclear polarisation transfer mechanisms such as the solid effect, cross effect, NOVEL, XiX, TOP, BEAM, and TPPM variants. The mathematics combines driven spin dynamics, relaxation, powder/MAS averaging, and steady-state or transient propagation.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Implementation structure

- E-15N cross-polarization experiment in the doubly rotating
- frame. Static powder simulation.
- Calculation time: seconds
- System specification
- Interactions
- Basis set
- Spinach housekeeping
- Experiment parameters
- Simulation
- Time axis generation
- Plotting

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `operator()`, `state()`, `powder()`, `cumsum()`, `kfigure()`, `kylabel()`, `kxlabel()`.
