# kernel/state.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/state.m`
- Signature: `rho=state(spin_system,states,spins,method)`
- Total lines: 321

## Purpose

Generates Hilbert space density matrices and Liouville space state vectors from their human-readable descriptions. Syntax: rho=state(spin_system,states,spins,method)

## Physical / mathematical content

- This file belongs to the `kernel` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content

- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 77-78: Default is to use consistent state norms; implemented by `if ~exist('method','var'), method='exact'; end`.
- Lines 80-81: In wavefunction space, empty set here; implemented by `if ~exist('spins','var'), spins=[]; end`.
- Lines 83-84: Check consistency; implemented by `grumble(spin_system,states,spins,method)`.
- Lines 86-87: Get the unit state; implemented by `switch spin_system.bas.formalism`.
- Lines 91-94: Unit population of T(0,0) state, normalisation is such because prod(spin_system.comp.mults) can be- come too large for double precision arithmetic; implemented by `unit=sparse(1,1,1,size(spin_system.bas.basis,1),1)`.
- Lines 98-100: Stretched unit matrix, normalisation matched to the Hilbert space because systems are small; implemented by `unit=speye(prod(spin_system.comp.mults)); unit=unit(:)`.
- Lines 104-105: Decide how to proceed; implemented by `switch spin_system.bas.formalism`.
- Lines 109-110: Choose the state vector generation methos; implemented by `switch method`.
- Lines 112-113: Careless normalisation; implemented by `case 'cheap'`.
- Lines 115-116: Parse the specification; implemented by `[opspecs,coeffs]=human2opspec(spin_system,states,spins)`.
- Lines 118-119: Compute correlation orders; implemented by `correlation_orders=sum(logical(spin_system.bas.basis),2)`.
- Lines 121-122: Locate each operator in the basis; implemented by `indices=zeros(size(coeffs))`.
- Lines 125-126: Find states with the same correlation order; implemented by `possibilities=(correlation_orders==nnz(opspecs{n}))`.
- Lines 128-129: Pin down the required state; implemented by `for k=find(opspecs{n})`.
- Lines 133-134: Double-check; implemented by `if nnz(possibilities)>1`.
- Lines 140-141: Locate the state; implemented by `indices(n)=find(possibilities)`.
- Lines 145-146: Assemble the state vector; implemented by `nrows=size(spin_system.bas.basis,1); ncols=1`.
- Lines 149-150: Careful normalisation; implemented by `case 'exact'`.

### Control flow inferred from the code

- Line 78: conditional branch on `~exist('method','var'), method='exact'; end`.
- Line 81: conditional branch on `~exist('spins','var'), spins=[]; end`.
- Line 87: dispatches on `spin_system.bas.formalism`; cases `'sphten-liouv'`, `'zeeman-liouv'`, `'sphten-liouv'`, `'cheap'`, `'exact'`, `'chem'`.
- Line 105: dispatches on `spin_system.bas.formalism`; cases `'sphten-liouv'`, `'cheap'`, `'exact'`, `'chem'`.
- Line 110: dispatches on `method`; cases `'cheap'`, `'exact'`, `'chem'`.
- Line 123: `parfor` loop over `n=1:numel(opspecs)`.
- Line 129: `for` loop over `k=find(opspecs{n})`.
- Line 134: conditional branch on `nnz(possibilities)>1`.
- Line 168: `for` loop over `n=1:numel(opspecs)`.
- Line 175: `for` loop over `k=1:numel(active_spins)`.
- Line 180: conditional branch on `nnz(species)~=1`.
- Line 219: `for` loop over `n=1:spin_system.comp.nspins`.

### Key state/data transformations

- Lines 94: computes `unit` using `unit=sparse(1,1,1,size(spin_system.bas.basis,1),1)`.
- Lines 116: computes `[opspecs,coeffs]` using `[opspecs,coeffs]=human2opspec(spin_system,states,spins)`.
- Lines 119: computes `correlation_orders` using `correlation_orders=sum(logical(spin_system.bas.basis),2)`.
- Lines 122: computes `indices` using `indices=zeros(size(coeffs))`.
- Lines 141: computes `indices(n)` using `indices(n)=find(possibilities)`.
- Lines 146: computes `nrows` using `nrows=size(spin_system.bas.basis,1); ncols=1`.
- Lines 147: computes `rho` using `rho=sparse(indices,ones(size(indices)),coeffs,nrows,ncols)`.
- Lines 165: computes `matrix_dim` using `matrix_dim=size(spin_system.bas.basis,1)`.
- Lines 171: computes `active_spins` using `active_spins=find(opspecs{n})`.
- Lines 174: computes `species` using `species=true(1,numel(spin_system.chem.parts))`.
- Lines 185: computes `coeffs(n)` using `coeffs(n)=coeffs(n)*spin_system.chem.concs(species)`.
- Lines 188: computes `A` using `A=superop(spin_system,opspecs{n},'left')`.
- Lines 216: computes `psi` using `psi=1`.
- Lines 222: computes `current_mult` using `current_mult=spin_system.comp.mults(n)`.
- Lines 223: computes `current_spin` using `current_spin=(current_mult-1)/2`.
- Lines 226: computes `levels` using `levels=fliplr((-current_spin):(current_spin))`.

### Local helper functions

- Line 247: `grumble()` — `function grumble(spin_system,states,spins,method)`. Check the formalism
  - Representative operation: `if (~isfield(spin_system,'bas'))||(~isfield(spin_system.bas,'formalism'))`.
  - Representative operation: `error('basis set information is missing, run basis() before calling this function.')`.

## Parameters / inputs

- 1. If states is a string and spins is a string
- states='Lz'; spins='13C';
- the function returns the sum of the corresponding single-spin densi-
- ty matrices (Hilbert space) or state vectors (Liouville space) on
- all spins of that type. Valid labels for states in this type of call
- are 'E' (identity), 'Lz', 'Lx', 'Ly', 'L+', 'L-', 'Tl,m' (irreduci-
- ble spherical tensor, l and m are integers), 'CTx', 'CTy', 'CTz',
- 'CT+','CT-' (central transition operators in the Zeeman basis). Va-
- lid labels for spins are standard isotope names, as well as 'elect-
- rons', 'nuclei', and 'all'.
- 2. If states is a string and spins is a vector
- states='Lz'; spins=[1 2 4];
- the function returns the sum of all single-spin density matrices
- (Hilbert space) or state vectors (Liouville space) for all spins
- with the specified numbers. Valid labels for states are the same as
- in Item 1 above.
- 3. If states is a cell array of strings and spins is a cell array
- of numbers:
- states={'Lz','L+'}; spins={1,2};
- then a product state density matrix (Hilbert space) or state vector
- (Liouville space) is produced. In the case above, Spinach will gene-
- rate LzS+ density matrix in Hilbert space or its state vector in Li-
- ouville space. Valid labels for operators are the same as in Item 1
- above.
- 4. For wavefunction formalism, states must be specified as an array
- of projection quantum numbers on all spins; in that case only two
- arguments are needed, for example, in a {'1H','1H','14N'} system:
- psi=state(spin_system,[-1/2 1/2 0])
- Method argument has the following effect in sphten-liouv formalism:
- 'cheap' -the state vector is generated without
- normalisation. For very large spin sys-
- tems this is much faster
- 'exact' -exact state vector with correct normalisation,
- this is the default when the last argument is
- skipped in the function call
- 'chem' -the exact state vector weighted with the
- concentrations specified in inter.chem.concs
- field under chemical kinetics parameters
- This option is ignored in zeeman-hilb and zeeman-liouv formalisms
- because there are no cheap shortcuts and kinetics is not available.

## Outputs

- rho -a Hilbert space density matrix or a Liouville
- space state vector

## Implementation structure

- Generates Hilbert space density matrices and Liouville space state
- vectors from their human-readable descriptions. Syntax:
- rho=state(spin_system,states,spins,method)
- 1. If states is a string and spins is a string
- states='Lz'; spins='13C';
- the function returns the sum of the corresponding single-spin densi-
- ty matrices (Hilbert space) or state vectors (Liouville space) on
- all spins of that type. Valid labels for states in this type of call
- are 'E' (identity), 'Lz', 'Lx', 'Ly', 'L+', 'L-', 'Tl,m' (irreduci-
- ble spherical tensor, l and m are integers), 'CTx', 'CTy', 'CTz',
- 'CT+','CT-' (central transition operators in the Zeeman basis). Va-
- lid labels for spins are standard isotope names, as well as 'elect-

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `exist()`, `grumble()`, `speye()`, `unit()`, `human2opspec()`, `logical()`, `nnz()`, `and()`, `indices()`, `operator()`, `spalloc()`, `true()`, `cellfun()`, `ismember()`, `active_spins()`, `coeffs()`.
