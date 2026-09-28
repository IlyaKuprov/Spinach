# kernel/states/partner_state.m

- Signature: `[A,descr]=partner_state(spin_system,set_spin,partners)`

## Purpose

Expands specified states of set_spin by Kronecker products with all combinations of the specified partner states. Syntax: `[A,descr]=partner_state(spin_system,set_spin,partners)`

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
