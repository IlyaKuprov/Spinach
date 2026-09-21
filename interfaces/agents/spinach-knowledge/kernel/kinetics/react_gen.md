# kernel/kinetics/react_gen.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/kinetics/react_gen.m`
- Signature: `G=react_gen(spin_system,reaction)`
- Total lines: 165

## Purpose

Chemical reaction generator builder. Syntax: G=react_gen(spin_system,reaction)

## Physical / mathematical content

- This file belongs to the `kernel` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content

- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

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

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `report()`, `cellfun()`, `all()`, `ismember()`, `drain_gen_idx()`, `destin_state()`, `source_state()`, `fill_gen_idx()`, `gen_idx()`, `complex()`, `num2str()`, `toc()`, `intersect()`, `cell2mat()`, `setdiff()`.
