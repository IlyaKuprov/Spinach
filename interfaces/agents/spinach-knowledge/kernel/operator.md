# kernel/operator.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/operator.m`
- Signature: `A=operator(spin_system,operators,spins,operator_type,format)`
- Total lines: 271

## Purpose

Generates Hilbert space operators or Liouville space superoperators from their human-readable descriptions. Syntax: A=operator(spin_system,operators,spins,operator_type,format)

## Physical / mathematical content

- This file belongs to the `kernel` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content

- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `numel()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 73-74: The default superoperator type is commutation; implemented by `if ~exist('operator_type','var'), operator_type='comm'; end`.
- Lines 76-77: The default sparse matrix format is CSC; implemented by `if ~exist('format','var'), format='csc'; end`.
- Lines 79-80: Check consistency; implemented by `grumble(spin_system,operators,spins,operator_type,format); tic`.
- Lines 82-83: Load the cache record if one exists; implemented by `if ismember('op_cache',spin_system.sys.enable)`.
- Lines 85-88: Combine specification, isotopes, and basis hash; implemented by `op_hash=md5_hash({operators,spins,operator_type, format,spin_system.comp.iso_hash, spin_system.bas.basis_hash})`.
- Lines 91-92: Get ValueStore; implemented by `if ~isworkernode`.
- Lines 98-99: Try to retrieve the operator from the ValueStore; implemented by `if isKey(store,op_hash), A=store(op_hash); return; end`.
- Lines 103-104: Parse the human specification into Spinach notation; implemented by `[opspecs,coeffs]=human2opspec(spin_system,operators,spins)`.
- Lines 106-107: Start a cell array; implemented by `A=cell(numel(opspecs),1)`.
- Lines 109-110: Decide how to proceed; implemented by `switch spin_system.bas.formalism`.
- Lines 112-113: Spherical tensor basis; implemented by `case 'sphten-liouv'`.
- Lines 115-116: Build summation terms; implemented by `parfor n=1:numel(opspecs)`.
- Lines 118-119: Get the superoperator; implemented by `A{n}=superop(spin_system,opspecs{n},operator_type)`.
- Lines 121-122: Apply the coefficient; implemented by `A{n}(:,3)=coeffs(n)*A{n}(:,3)`.
- Lines 126-129: Zeeman basis formalisms; implemented by `case {'zeeman-wavef', 'zeeman-hilb', 'zeeman-liouv'}`.
- Lines 131-132: Parallelisation efficiency; implemented by `mults=spin_system.comp.mults`.
- Lines 138-139: Start the kron; implemented by `B=sparse(coeffs(n))`.
- Lines 141-142: Over particles; implemented by `for k=1:numel(opspecs{n})`.

### Control flow inferred from the code

- Line 74: conditional branch on `~exist('operator_type','var'), operator_type='comm'; end`.
- Line 77: conditional branch on `~exist('format','var'), format='csc'; end`.
- Line 83: conditional branch on `ismember('op_cache',spin_system.sys.enable)`.
- Line 92: conditional branch on `~isworkernode`.
- Line 99: conditional branch on `isKey(store,op_hash), A=store(op_hash); return; end`.
- Line 110: dispatches on `spin_system.bas.formalism`; cases `'sphten-liouv'`, `{'zeeman-wavef',`.
- Line 116: `parfor` loop over `n=1:numel(opspecs)`.
- Line 136: `parfor` loop over `n=1:numel(opspecs)`.
- Line 142: `for` loop over `k=1:numel(opspecs{n})`.
- Line 156: conditional branch on `strcmp(formalism,'zeeman-liouv')`.
- Line 176: conditional branch on `strcmp(format,'csc')`.
- Line 179: dispatches on `spin_system.bas.formalism`; cases `'sphten-liouv'`, `{'zeeman-wavef','zeeman-hilb'}`, `'zeeman-liouv'`.
- Line 209: conditional branch on `ismember('op_cache',spin_system.sys.enable)`.
- Line 212: conditional branch on `~isKey(store,op_hash)`.

### Key state/data transformations

- Lines 86-88: computes `op_hash` using `op_hash=md5_hash({operators,spins,operator_type, format,spin_system.comp.iso_hash, spin_system.bas.basis_hash})`.
- Lines 93: computes `store` using `store=gcp('nocreate').ValueStore`.
- Lines 104: computes `[opspecs,coeffs]` using `[opspecs,coeffs]=human2opspec(spin_system,operators,spins)`.
- Lines 107: computes `A` using `A=cell(numel(opspecs),1)`.
- Lines 119: computes `A{n}` using `A{n}=superop(spin_system,opspecs{n},operator_type)`.
- Lines 122: computes `A{n}(:,3)` using `A{n}(:,3)=coeffs(n)*A{n}(:,3)`.
- Lines 132: computes `mults` using `mults=spin_system.comp.mults`.
- Lines 133: computes `formalism` using `formalism=spin_system.bas.formalism`.
- Lines 139: computes `B` using `B=sparse(coeffs(n))`.
- Lines 145: computes `[L,M]` using `[L,M]=lin2lm(opspecs{n}(k))`.
- Lines 148: computes `IST` using `IST=irr_sph_ten(mults(k),L)`.
- Lines 161: computes `[rows,cols,vals]` using `[rows,cols,vals]=find(B); A{n}=[rows cols vals]`.
- Lines 184: computes `matrix_dim` using `matrix_dim=size(spin_system.bas.basis,1)`.

### Local helper functions

- Line 221: `grumble()` — `function grumble(spin_system,operators,spins,operator_type,format)`.
  - Representative operation: `if ~isfield(spin_system,'bas')`.
  - Representative operation: `error('basis set information is missing, run basis() before calling this function.')`.

## Parameters / inputs

- 1. If operators is a string and spins is a string
- operators='Lz'; spins='13C';
- the function returns the sum of the corresponding single-spin operators
- (Hilbert space) or superoperators (Liouville space) on all spins of that
- type. Valid labels for states in this type of call are 'E' (identity),
- 'Lz', 'Lx', 'Ly', 'L+', 'L-', 'Tl,m' (irreducible spherical tensor, l
- and m are integers), 'CTx', 'CTy', 'CTz', 'CT+', 'CT-' (central transi-
- tion operators in the Zeeman basis). Valid labels for spins are standard
- isotope names, as well as 'electrons', 'nuclei', and 'all'.
- 2. If operators is a string and spins is a vector
- operators='Lz'; spins=[1 2 4];
- the function returns the sum of all single-spin operators (Hilbert space)
- or superoperators (Liouville space) for all spins with the specified num-
- bers. Valid labels for operators are the same as in Item 1 above.
- 3. If operators is a cell array of strings and spins is a cell array of
- numbers
- operators={'Lz','L+'}; spins={1,2};
- then a product operator (Hilbert space) or its superoperator (Liouville
- space) is produced. In the case above, Spinach will generate LzS+ in Hil-
- bert space or its specified superoperator in Liouville space. Valid la-
- bels for operators are the same as in Item 1 above.
- In Liouville space calculations, operator_type can be set to:
- 'left' -produces left side product superoperator
- 'right' -produces right side product superoperator
- 'comm' -produces commutation superoperator (default)
- 'acomm' -produces anticommutation superoperator
- In Hilbert space calculations operator_type parameter is ignored, and the
- operator itself is always returned.
- The format parameter refers to the format of the output: 'csc' returns a
- Matlab sparse matrix, 'xyz' returns a [rows, cols, vals] array.

## Outputs

- A -a CSC sparse (default) or a [rows, cols, vals] repre-
- sentation of a spin operator or superoperator.
- Notes: WARNING -a product of two commutation superoperators is NOT a com-
- mutation superoperator of a product. In Liouville space, you cannot
- generate single-spin superoperators and multiply them up.
- Note: operator caching is supported, add 'op_cache' to sys.enable array
- to enable; make sure your scratch storage is fast.

## Implementation structure

- Generates Hilbert space operators or Liouville space superoperators from
- their human-readable descriptions. Syntax:
- A=operator(spin_system,operators,spins,operator_type,format)
- 1. If operators is a string and spins is a string
- operators='Lz'; spins='13C';
- the function returns the sum of the corresponding single-spin operators
- (Hilbert space) or superoperators (Liouville space) on all spins of that
- type. Valid labels for states in this type of call are 'E' (identity),
- 'Lz', 'Lx', 'Ly', 'L+', 'L-', 'Tl,m' (irreducible spherical tensor, l
- and m are integers), 'CTx', 'CTy', 'CTz', 'CT+', 'CT-' (central transi-
- tion operators in the Zeeman basis). Valid labels for spins are standard
- isotope names, as well as 'electrons', 'nuclei', and 'all'.

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `exist()`, `grumble()`, `ismember()`, `md5_hash()`, `gcp()`, `getCurrentValueStore()`, `isKey()`, `store()`, `human2opspec()`, `superop()`, `coeffs()`, `lin2lm()`, `irr_sph_ten()`, `mults()`, `strcmp()`, `hilb2liouv()`.
