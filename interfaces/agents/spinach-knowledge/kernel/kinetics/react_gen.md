# kernel/kinetics/react_gen.m

- Signature: `G=react_gen(spin_system,reaction)`

## Purpose

Builds state-space generators for the specified chemical reaction, mapping reactant basis states to corresponding product states.

## Physical / mathematical content

For each reactant, the generator removes the source state and adds the matched product state when that product state is present in the basis. The mapping is defined by `reaction.matching`.

## Numerical / algorithmic content

The routine scans basis states, identifies their chemical part, and constructs sparse complex matrices. It checks that reactants and products are disjoint, have the same total number of spins, and agree with the matching map.

## Parameters / inputs

- `spin_system` - Spinach system description, including the chemical parts and basis.
- `reaction.reactants` - vector of indices for reactants among `spin_system.chem.parts`.
- `reaction.products` - vector of indices for products among `spin_system.chem.parts`.
- `reaction.matching` - two-column matrix mapping spin indices from the reactant side (left column) to the product side (right column).

## Outputs

- `G` - cell array with one sparse generator matrix per reactant, representing source-state drainage and matching product-state filling in the basis.

## Implementation structure

The routine validates the reaction specification, scans the basis in parallel, records drainage for each reactant state and filling for an available matched product state, then assembles one sparse matrix per reactant.
