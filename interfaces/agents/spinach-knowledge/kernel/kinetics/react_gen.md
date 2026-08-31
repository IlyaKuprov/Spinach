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

## Code-derived implementation details

### Comment-guided execution stages

- Lines 33-34: Check consistency; implemented by `grumble(spin_system,reaction)`.
- Lines 36-37: Inform the user and get the timer going; implemented by `report(spin_system,'building reaction generators ')`.
- Lines 40-42: Preallocate reaction generator arrays to be [reactant destin source coeff]; implemented by `nstates=size(spin_system.bas.basis,1)`.
- Lines 46-47: Loop over the basis set; implemented by `parfor n=1:nstates`.
- Lines 49-50: Extract the state; implemented by `source_state=spin_system.bas.basis(n,:)`.
- Lines 52-53: Find participating spins; implemented by `[~,spins_involved]=find(source_state)`.
- Lines 55-56: Build reaction generators; implemented by `if ~isempty(spins_involved)`.
- Lines 58-60: Determine the host substance; implemented by `host_subst=cellfun(@(x)all(ismember(spins_involved,x)), spin_system.chem.parts)`.
- Lines 63-64: Double-check basis state indexing; implemented by `if numel(host_subst)~=1, error('basis set indexing problem.'); end`.
- Lines 66-67: Determine host substance type; implemented by `this_is_reactants=ismember(host_subst,reaction.reactants)`.
- Lines 70-71: Double-check indexing; implemented by `if this_is_reactants&&this_is_products`.
- Lines 75-76: Only reactants react; implemented by `if this_is_reactants`.
- Lines 78-79: Add to reactant drain generator; implemented by `idx=find(reaction.reactants==host_subst)`.
- Lines 82-83: Build the destination state; implemented by `destin_state=zeros(1,numel(source_state))`.
- Lines 87-88: Look for the destination state in the basis set and double-check indexing; implemented by `[destin_exists,destin_index]=ismember(destin_state,spin_system.bas.basis,'rows')`.
- Lines 91-92: Build product fill generator; implemented by `if destin_exists`.
- Lines 94-95: Find participating spins; implemented by `[~,spins_involved]=find(destin_state)`.
- Lines 102-103: Double-check indexing; implemented by `if numel(host_subst)~=1`.

### Control flow inferred from the code

- Line 47: `parfor` loop over `n=1:nstates`.
- Line 56: conditional branch on `~isempty(spins_involved)`.
- Line 64: conditional branch on `numel(host_subst)~=1, error('basis set indexing problem.'); end`.
- Line 71: conditional branch on `this_is_reactants&&this_is_products`.
- Line 76: conditional branch on `this_is_reactants`.
- Line 89: conditional branch on `numel(destin_index)>1, error('invalid basis set specification'); end`.
- Line 92: conditional branch on `destin_exists`.
- Line 103: conditional branch on `numel(host_subst)~=1`.
- Line 124: `for` loop over `n=1:numel(reaction.reactants)`.

### Key state/data transformations

- Lines 38: computes `timer_react_gen` using `timer_react_gen=tic`.
- Lines 42: computes `nstates` using `nstates=size(spin_system.bas.basis,1)`.
- Lines 43: computes `drain_gen_idx` using `drain_gen_idx=zeros(nstates,4)`.
- Lines 44: computes `fill_gen_idx` using `fill_gen_idx=zeros(nstates,4)`.
- Lines 50: computes `source_state` using `source_state=spin_system.bas.basis(n,:)`.
- Lines 53: computes `[~,spins_involved]` using `[~,spins_involved]=find(source_state)`.
- Lines 59-60: computes `host_subst` using `host_subst=cellfun(@(x)all(ismember(spins_involved,x)), spin_system.chem.parts)`.
- Lines 61: computes `[~,host_subst]` using `[~,host_subst]=find(host_subst)`.
- Lines 67: computes `this_is_reactants` using `this_is_reactants=ismember(host_subst,reaction.reactants)`.
- Lines 68: computes `this_is_products` using `this_is_products=ismember(host_subst,reaction.products)`.
- Lines 80: computes `drain_gen_idx(n,:)` using `drain_gen_idx(n,:)=[idx n n -1]`.
- Lines 83: computes `destin_state` using `destin_state=zeros(1,numel(source_state))`.
- Lines 84: computes `destin_state(reaction.matching(:,2))` using `destin_state(reaction.matching(:,2))=source_state(reaction.matching(:,1))`.
- Lines 88: computes `[destin_exists,destin_index]` using `[destin_exists,destin_index]=ismember(destin_state,spin_system.bas.basis,'rows')`.
- Lines 108: computes `fill_gen_idx(n,:)` using `fill_gen_idx(n,:)=[idx destin_index n 1]`.
- Lines 119: computes `gen_idx` using `gen_idx=[drain_gen_idx; fill_gen_idx]`.
- Lines 123: computes `G` using `G=cell([numel(reaction.reactants) 1])`.
- Lines 128: computes `G{n}` using `G{n}=complex(G{n})`.

### Local helper functions

- Line 138: `grumble()` — `function grumble(spin_system,reaction)`.
  - Representative operation: `if ~isempty(intersect(reaction.reactants, reaction.products))`.
  - Representative operation: `reaction.products))`.

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
