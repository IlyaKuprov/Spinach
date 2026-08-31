# kernel/utilities/human2opspec.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/utilities/human2opspec.m`
- Signature: `[opspecs,coeffs]=human2opspec(spin_system,operators,spins)`
- Total lines: 374

## Purpose

Converts user-friendly descriptions of spin states and operators into the formal description (opspec) used by Spinach kernel. Syntax: [opspecs,coeffs]=human2opspec(spin_system,operators,spins)

## Physical / mathematical content

- General mathematical and infrastructure utilities. This area contains finite differences, perturbation theory, graph algorithms, spectral densities, tensor algebra, hash/report helpers, and other reusable numerical components.
- The effective hardware model is a weakly anharmonic oscillator. Duffing nonlinearity breaks equal level spacing and allows qubit-like addressability within a truncated bosonic ladder.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `numel()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 57-58: Check input validity; implemented by `grumble(operators,spins)`.
- Lines 60-61: 'Lz','1H' type call returns a sum; implemented by `if ischar(operators)&&ischar(spins)`.
- Lines 63-64: Parse the specification; implemented by `switch spins`.
- Lines 68-69: Use all spins in the system; implemented by `spin_numbers=1:spin_system.comp.nspins`.
- Lines 73-74: Include electrons of any multiplicity; implemented by `spin_numbers=find(cellfun(@(x)strncmp(x,'E',1),spin_system.comp.isotopes))`.
- Lines 78-79: Include non-electrons; implemented by `spin_numbers=find(~cellfun(@(x)strncmp(x,'E',1),spin_system.comp.isotopes))`.
- Lines 83-84: Count specific spins; implemented by `spin_numbers=find(strcmp(spin_system.comp.isotopes,spins))`.
- Lines 88-89: Bomb out if the list of spins is empty; implemented by `if numel(spin_numbers)==0, error('no such spins in the system.'); end`.
- Lines 91-92: Preallocate descriptor arrays; implemented by `coeffs=cell(numel(spin_numbers),1); opspecs=cell(numel(spin_numbers),1)`.
- Lines 94-95: Run recursive calls; implemented by `for n=1:numel(spin_numbers)`.
- Lines 99-100: Concatenate descriptor arrays and return the outcome; implemented by `opspecs=vertcat(opspecs{:}); coeffs=cell2mat(coeffs); return`.
- Lines 102-103: 'Lz',[1 2 4] type call returns a sum; implemented by `elseif ischar(operators)&&isnumeric(spins)`.
- Lines 105-106: Preallocate descriptor arrays; implemented by `coeffs=cell(numel(spins),1); opspecs=cell(numel(spins),1)`.
- Lines 108-109: Run recursive calls; implemented by `for n=1:numel(spins)`.
- Lines 116-117: {'Lz','L+'},{1,4} type call returns an operator; implemented by `elseif iscell(operators)&&iscell(spins)`.
- Lines 119-120: Start with an empty opspec and a unit coefficient; implemented by `opspecs=zeros(1,spin_system.comp.nspins); coeffs=1`.
- Lines 122-123: Parse operator type; implemented by `for n=1:numel(operators)`.
- Lines 125-126: Different for spins and bosons; implemented by `switch spin_system.comp.types{spins{n}}`.

### Control flow inferred from the code

- Line 61: conditional branch on `ischar(operators)&&ischar(spins)`.
- Line 64: dispatches on `spins`; cases `'all'`, `'electrons'`, `'nuclei'`.
- Line 89: conditional branch on `numel(spin_numbers)==0, error('no such spins in the system.'); end`.
- Line 95: `for` loop over `n=1:numel(spin_numbers)`.
- Line 109: `for` loop over `n=1:numel(spins)`.
- Line 123: `for` loop over `n=1:numel(operators)`.
- Line 126: dispatches on `spin_system.comp.types{spins{n}}`; cases `'S'`, `'E'`, `'L+'`, `'L-'`, `'Lx'`, `'Ly'`, `'Lz'`, `'CTx'`, `'CTy'`, `'CTz'`, ….
- Line 131: dispatches on `operators{n}`; cases `'E'`, `'L+'`, `'L-'`, `'Lx'`, `'Ly'`, `'Lz'`, `'CTx'`, `'CTy'`, `'CTz'`, `'CT+'`, ….
- Line 218: conditional branch on `strcmp(operators{n}(1),'T')`.
- Line 221: conditional branch on `isempty(regexp(operators{n},'^T([\+\-]?\d+),([\+\-]?\d+)$','once'))`.
- Line 229: conditional branch on `(l<0)||(abs(m)>l)||(mod(l,1)~=0)||(mod(m,1)~=0)||`.
- Line 241: conditional branch on `isempty(regexp(operators{n},'^ZL([1-9]\d*)$','once'))`.
- Line 250: conditional branch on `(level_number<1)||(mod(level_number,1)~=0)||`.
- Line 276: conditional branch on `strcmp(operators{n},'E')`.

### Key state/data transformations

- Lines 69: computes `spin_numbers` using `spin_numbers=1:spin_system.comp.nspins`.
- Lines 92: computes `coeffs` using `coeffs=cell(numel(spin_numbers),1); opspecs=cell(numel(spin_numbers),1)`.
- Lines 96: computes `[opspecs{n},coeffs{n}]` using `[opspecs{n},coeffs{n}]=human2opspec(spin_system,{operators},{spin_numbers(n)})`.
- Lines 100: computes `opspecs` using `opspecs=vertcat(opspecs{:}); coeffs=cell2mat(coeffs); return`.
- Lines 136: computes `opspecs(:,spins{n})` using `opspecs(:,spins{n})=0`.
- Lines 157: computes `opspecs_a` using `opspecs_a=opspecs; opspecs_a(:,spins{n})=1`.
- Lines 158: computes `opspecs_b` using `opspecs_b=opspecs; opspecs_b(:,spins{n})=3`.
- Lines 178: computes `[ct_states,ct_coeffs]` using `[ct_states,ct_coeffs]=ct2ist(spin_system.comp.mults(spins{n}),'x')`.
- Lines 179: computes `states` using `states=kron(ct_states,ones(size(opspecs,1),1))`.
- Lines 226: computes `indices` using `indices=textscan(operators{n},'T%n,%n'); l=indices{1}; m=indices{2}`.
- Lines 246: computes `level_number` using `level_number=textscan(operators{n},'ZL%n')`.
- Lines 256: computes `[zl_states,zl_coeffs]` using `[zl_states,zl_coeffs]=enlev2ist(spin_system.comp.mults(spins{n}),level_number,'S')`.
- Lines 273: computes `nlevels` using `nlevels=spin_system.comp.mults(spins{n})`.
- Lines 298: computes `[b_states,b_coeffs]` using `[b_states,b_coeffs]=enlev2ist(spin_system.comp.mults(spins{n}),level_number,'B')`.

### Local helper functions

- Line 335: `grumble()` — `function grumble(operators,spins)`.
  - Representative operation: `if (~(ischar(operators)&&ischar(spins)))&& (~(iscell(operators)&&iscell(spins)))&& (~(ischar(operators)&&isnumeric(spins)))`.
  - Representative operation: `(~(iscell(operators)&&iscell(spins)))&& (~(ischar(operators)&&isnumeric(spins)))`.

## Parameters / inputs

- 1. If both inputs are strings
- [opspecs,coeffs]=human2opspec(spin_system,'Lz','13C')
- the function returns a list of single-spin opspecs for all spins with the
- specified name. In the example above, the list of Lz operator specificati-
- ons for all 13C nuclei in the system would be returned. Valid labels for
- states in this type of call are 'E' (identity), 'Lz', 'Lx', 'Ly', 'L+',
- 'L-', 'Tl,m' (irreducible spherical tensor, l and m are integers), 'CTx',
- 'CTy', 'CTz', 'CT+','CT-' (central transition operators in the Zeeman ba-
- sis). Valid labels for spins are standard isotope names, as well as 'elec-
- trons', 'nuclei', and 'all'.
- 2. If one input is a string and the other is a vector
- [opspecs,coeffs]=human2opspec(spin_system,'Lz',[1 2 4])
- the function returns a list of single-spin opspecs for all spins with the
- specified number. In the example above, the list of Lz operator specifica-
- tions for all 13C nuclei in the system would be returned. Valid labels for
- operators are the same as in Item 1 above.
- 3. If the two inputs are a cell array of strings and a cell array of num-
- bers, a product operator specification is produced
- [opspecs,coeffs]=human2opspec(spin_system,{'Lz','Ly'},{1,2})
- would return the Lz(x)Ly product operator specification with Lz on spin 1
- and Ly on spin 2. Valid labels for operators are the same as in Item 1.

## Outputs

- opspecs -Spinach operator specification: a cell array of
- row vectors specifying which operator enters the
- Kronecker product for which spin.
- coeffs -coefficient with which each of the Kronecker pro-
- ducts enters the overall sum.
- Notes: direct calls to this function are not necessary, use operator.m and
- state.m functions instead.

## Implementation structure

- Converts user-friendly descriptions of spin states and operators into the
- formal description (opspec) used by Spinach kernel. Syntax:
- [opspecs,coeffs]=human2opspec(spin_system,operators,spins)
- 1. If both inputs are strings
- [opspecs,coeffs]=human2opspec(spin_system,'Lz','13C')
- the function returns a list of single-spin opspecs for all spins with the
- specified name. In the example above, the list of Lz operator specificati-
- ons for all 13C nuclei in the system would be returned. Valid labels for
- states in this type of call are 'E' (identity), 'Lz', 'Lx', 'Ly', 'L+',
- 'L-', 'Tl,m' (irreducible spherical tensor, l and m are integers), 'CTx',
- 'CTy', 'CTz', 'CT+','CT-' (central transition operators in the Zeeman ba-
- sis). Valid labels for spins are standard isotope names, as well as 'elec-

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `ischar()`, `cellfun()`, `strncmp()`, `strcmp()`, `spin_numbers()`, `vertcat()`, `cell2mat()`, `spins()`, `iscell()`, `opspecs()`, `opspecs_a()`, `opspecs_b()`, `ct2ist()`, `regexp()`, `textscan()`.
