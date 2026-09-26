# kernel/kinetics/react_gen.m

- Signature: `G=react_gen(spin_system,reaction)`

## Purpose

Chemical reaction generator builder. Syntax: G=react_gen(spin_system,reaction)

## Physical / mathematical content

## Numerical / algorithmic content

- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.

## Parameters / inputs

- reaction.reactants -a vector of integers specifying
- which parts declared in the in-
- put (chem.parts) are reactants
- reaction.products -a vector of integers specifying
- which parts declared in the in-
- put (chem.parts) are products
- reaction.matching -a matrix with two columns, spe-
- cifying which spin in the reac-
- tants list (left column) becom-
- es which spin in the product
- list (right column)

## Outputs

- G -a cell array of matrices, one per reactant, map-
- ping each state of the reactant state space into
- its destination in the product state space

## Implementation structure

- Chemical reaction generator builder. Syntax:
- G=react_gen(spin_system,reaction)
- reaction.reactants -a vector of integers specifying
- which parts declared in the in-
- put (chem.parts) are reactants
- reaction.products -a vector of integers specifying
- put (chem.parts) are products
- reaction.matching -a matrix with two columns, spe-
- cifying which spin in the reac-
- tants list (left column) becom-
- es which spin in the product
- list (right column)
