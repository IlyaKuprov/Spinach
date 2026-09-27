# examples/microfluidics/reacting_flow.m

- Signature: `reacting_flow()`

## Purpose

Flow in the absence of spin dynamics, but presence of two unidirectional second-order chemical reactions. Simulation time: seconds.

## Physical / mathematical content

- This example has no spin dynamics: it combines flow and strong diffusion on an imported COMSOL mesh with two competing second-order cycloaddition reactions. The local chemistry is evaluated from reactant concentrations in each mesh cell.

## Numerical / algorithmic content

- At each time step, the code builds a cell-specific chemical generator, combines it with the mesh flow/diffusion generator, and advances the flattened concentration trajectory with `step`.
- A `parfor` loop evaluates the local reaction generator independently for each mesh cell; the four plotted fields are the two reactants and two products.

## Implementation structure

- Flow in the absence of spin dynamics, but presence of two
- unidirectional second-order chemical reactions.
- Simulation time: seconds.
- Import hydrodynamics information
- No spin system here
- Rate constants, mol/(L*s)
- Cycloaddition reaction generator, including solvent
- Strong diffusion
- Timing parameters
- Get diffusion and flow generator
- Trajectory preallocation and the initial state
- Time evolution loop
