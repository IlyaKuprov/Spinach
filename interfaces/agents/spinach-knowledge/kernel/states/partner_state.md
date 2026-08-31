# kernel/states/partner_state.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/states/partner_state.m`
- Signature: `[A,descr]=partner_state(spin_system,set_spin,partners)`
- Total lines: 208

## Purpose

Partner state expansion; a given state of the specified spins is kroneckered with all combinations of the specified states of specfied partner spins. Syntax: [A,descr]=partner_state(spin_system,active_spin,partners)

## Physical / mathematical content

- State-construction utilities. These routines build equilibrium states, singlets, triplets, partner-state expansions, and physically meaningful density operators in the active basis.

## Numerical / algorithmic content

- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 71-72: Check consistency; implemented by `grumble(spin_system,set_spin,partners)`.
- Lines 74-75: Write the states of the spins that stay unchanged; implemented by `base_state=repmat({'E'},1,spin_system.comp.nspins)`.
- Lines 80-81: Pull out partner spin lists and partner state sets; implemented by `partner_spin_lists=cellfun(@(spec)spec{2},partners,'UniformOutput',false)`.
- Lines 84-85: Compile the active partner list; implemented by `active_partners=[partner_spin_lists{:}]`.
- Lines 87-90: Compile state lists for each partner; implemented by `partner_state_sets=cellfun(@(x,y)repmat({x},1,numel(y)), partner_state_sets, partner_spin_lists,'UniformOutput',false)`.
- Lines 93-94: Start descriptor; implemented by `descr={base_state}`.
- Lines 96-97: Over active partner spins; implemented by `for n=1:numel(active_partners)`.
- Lines 99-100: Pull current partner; implemented by `partner=active_partners(n)`.
- Lines 102-103: Over the partner state set; implemented by `for k=1:numel(partner_state_sets{n})`.
- Lines 105-106: Pull current partner spin state; implemented by `partner_state=partner_state_sets{n}{k}`.
- Lines 108-109: Scan descriptors; implemented by `for m=1:numel(descr)`.
- Lines 111-112: Pull descriptor line; implemented by `current_line=descr{m}`.
- Lines 114-115: Modify the descriptor line if necessary; implemented by `if ~strcmp(current_line{partner},partner_state)`.
- Lines 117-118: Write in the new state and append; implemented by `current_line{partner}=partner_state`.
- Lines 129-130: Get the full spin list; implemented by `full_spin_list=num2cell(1:spin_system.comp.nspins)`.
- Lines 132-133: Generate the states; implemented by `A=cell(1,numel(descr))`.

### Control flow inferred from the code

- Line 76: `for` loop over `n=1:numel(set_spin)`.
- Line 97: `for` loop over `n=1:numel(active_partners)`.
- Line 103: `for` loop over `k=1:numel(partner_state_sets{n})`.
- Line 109: `for` loop over `m=1:numel(descr)`.
- Line 115: conditional branch on `~strcmp(current_line{partner},partner_state)`.
- Line 134: `parfor` loop over `n=1:numel(descr)`.

### Key state/data transformations

- Lines 75: computes `base_state` using `base_state=repmat({'E'},1,spin_system.comp.nspins)`.
- Lines 77: computes `base_state{set_spin{n}{2}}` using `base_state{set_spin{n}{2}}=set_spin{n}{1}`.
- Lines 81: computes `partner_spin_lists` using `partner_spin_lists=cellfun(@(spec)spec{2},partners,'UniformOutput',false)`.
- Lines 82: computes `partner_state_sets` using `partner_state_sets=cellfun(@(spec)spec{1},partners,'UniformOutput',false)`.
- Lines 85: computes `active_partners` using `active_partners=[partner_spin_lists{:}]`.
- Lines 94: computes `descr` using `descr={base_state}`.
- Lines 100: computes `partner` using `partner=active_partners(n)`.
- Lines 106: computes `partner_state` using `partner_state=partner_state_sets{n}{k}`.
- Lines 112: computes `current_line` using `current_line=descr{m}`.
- Lines 118: computes `current_line{partner}` using `current_line{partner}=partner_state`.
- Lines 130: computes `full_spin_list` using `full_spin_list=num2cell(1:spin_system.comp.nspins)`.
- Lines 133: computes `A` using `A=cell(1,numel(descr))`.
- Lines 135: computes `A{n}` using `A{n}=state(spin_system,descr{n},full_spin_list)`.

### Local helper functions

- Line 141: `grumble()` — `function grumble(spin_system,set_spin,partners)`.
  - Representative operation: `if (~isfield(spin_system,'comp'))||(~isfield(spin_system.comp,'nspins'))|| (~isnumeric(spin_system.comp.nspins))||(~isreal(spin_system.comp.nspins))|| (~isscalar(spin_sy…`.
  - Representative operation: `(~isnumeric(spin_system.comp.nspins))||(~isreal(spin_system.comp.nspins))|| (~isscalar(spin_system.comp.nspins))||(spin_system.comp.nspins<1)|| (mod(spin_system.comp.nsp…`.

## Parameters / inputs

- set_spin -a cell array of two-element cell arrays
- with the first element giving the state
- of the spin and the second element num-
- ber on the isotope list, for example
- {{'L+',3}}
- These spins will have their states set
- immutably as specified.
- partners -a cell array of partner state specifica-
- tions of the form
- { {states_a,spins_a},...
- {states_b,spins_b},... }
- where states_a,b are cell arrays of sta-
- tes that the partner spins can have and
- spins_a,b are lists of numbers of those
- spins on the main isotope list, e.g.
- { {{'E' 'Lz'},[1 5]},...
- {{'L+' 'L-'},[2 7]},... }
- These spins will have their state varied
- combinatorially as specified.

## Outputs

- A -a cell array of spin states (matrices in Hilbert space,
- vectors in Liouville space) with the active spin in the
- specified state and the partner spins in all combinati-
- ons specified by the user. All spins not explicitly
- mentioned in the input will be in their 'E' states.
- descr -a cell array of product structure descriptors for
- each element of A
- Example: in a five-spin system, the following call
- A=partner_state(spin_system,{{'L+',2}},{{{'E','Lz'},[1 3]}})
- will return the following state array
- A={state(spin_system,{'E' ,'L+','E' ,'E' ,'E'},{1 2 3 4 5}),...
- state(spin_system,{'Lz','L+','E' ,'E' ,'E'},{1 2 3 4 5}),...
- state(spin_system,{'E' ,'L+','Lz','E' ,'E'},{1 2 3 4 5}),...
- state(spin_system,{'Lz','L+','Lz','E' ,'E'},{1 2 3 4 5})};
- and the following descriptor array
- descr={{'E' ,'L+','E' ,'E' ,'E'},...
- {'Lz','L+','E' ,'E' ,'E'},...
- {'E' ,'L+','Lz','E' ,'E'},...
- {'Lz','L+','Lz','E' ,'E'}};

## Implementation structure

- Partner state expansion; a given state of the specified spins
- is kroneckered with all combinations of the specified states
- of specfied partner spins. Syntax:
- [A,descr]=partner_state(spin_system,active_spin,partners)
- set_spin -a cell array of two-element cell arrays
- with the first element giving the state
- of the spin and the second element num-
- ber on the isotope list, for example
- {{'L+',3}}
- These spins will have their states set
- immutably as specified.
- partners -a cell array of partner state specifica-

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `cellfun()`, `active_partners()`, `strcmp()`, `num2cell()`, `state()`, `isfield()`, `isscalar()`, `iscell()`, `ischar()`, `set_spin_idx()`, `any()`, `isvector()`, `ismember()`.
