# kernel/operators/superop.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/operators/superop.m`
- Signature: `A=superop(spin_system,opspec,side)`
- Total lines: 212

## Purpose

Sided product superoperator in the spherical tensor basis set. Returns superoperators corresponding to right or left multiplication of a den- sity matrix by a user-specified operator. Syntax: A=superop(spin_system,opspec,side)

## Physical / mathematical content

- Operator-construction utilities. They build bases and irreducible tensor representations for spin, bosonic, and transition operators.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 38-39: Issue a recursive call if appropriate; implemented by `if strcmp(side,'comm')`.
- Lines 48-49: Validate the input; implemented by `grumble(spin_system,opspec)`.
- Lines 51-52: Determine the relevant spins; implemented by `active_spins=find(opspec)`.
- Lines 54-55: For unit operator use a shortcut; implemented by `if isempty(active_spins)`.
- Lines 62-63: Preallocate source state index; implemented by `source=cell(1,numel(active_spins))`.
- Lines 65-66: Preallocate destination state index; implemented by `destin=cell(1,numel(active_spins))`.
- Lines 68-69: Preallocate structure coefficients table; implemented by `struct=cell(1,numel(active_spins))`.
- Lines 71-72: Loop over the relevant spins; implemented by `for n=1:length(active_spins)`.
- Lines 74-75: Current spin multiplicity and state index; implemented by `mult=spin_system.comp.mults(active_spins(n))`.
- Lines 78-79: Extract pages corresponding to the current state; implemented by `switch side`.
- Lines 83-84: Extract left product pages from Lie structure tables; implemented by `pt=squeeze(spin_system.bas.lpst{mult}(table_idx,:,:))`.
- Lines 88-89: Extract right product pages from Lie structure tables; implemented by `pt=squeeze(spin_system.bas.rpst{mult}(table_idx,:,:))`.
- Lines 93-94: Complain and bomb out; implemented by `error('invalid side specification.')`.
- Lines 98-99: Convert product action table to indices; implemented by `[destin{n},source{n},struct{n}]=find(pt)`.
- Lines 101-102: Spinach uses 0 index for the unit matrix; implemented by `source{n}=source{n}-1; destin{n}=destin{n}-1`.
- Lines 106-107: Compute the structure coefficients for the relevant sub-algebra; implemented by `from=source{1}; to=destin{1}; coeff=struct{1}`.
- Lines 116-117: Lift the basis columns corresponding to the relevant spins; implemented by `basis_cols=spin_system.bas.basis(:,active_spins)`.
- Lines 119-120: For commutation superoperators remove commuting paths; implemented by `if ismember(side,{'leftofcomm','rightofcomm'})`.

### Control flow inferred from the code

- Line 39: conditional branch on `strcmp(side,'comm')`.
- Line 55: conditional branch on `isempty(active_spins)`.
- Line 72: `for` loop over `n=1:length(active_spins)`.
- Line 79: dispatches on `side`; cases `{'left','leftofcomm'}`, `{'right','rightofcomm'}`.
- Line 108: `for` loop over `n=2:numel(active_spins)`.
- Line 120: conditional branch on `ismember(side,{'leftofcomm','rightofcomm'})`.
- Line 129: `for` loop over `n=1:size(from,1)`.
- Line 133: `for` loop over `m=1:size(from,2)`.
- Line 144: conditional branch on `subsp_dim>0`.
- Line 148: `for` loop over `m=1:size(to,2)`.
- Line 156: conditional branch on `isequal(source_subsp,destin_subsp)`.
- Line 179: conditional branch on `isempty(A)`.

### Key state/data transformations

- Lines 40: computes `A` using `A=superop(spin_system,opspec,'leftofcomm')`.
- Lines 41: computes `B` using `B=superop(spin_system,opspec,'rightofcomm')`.
- Lines 42: computes `B(:,3)` using `B(:,3)=-B(:,3); A=[A; B]; return`.
- Lines 52: computes `active_spins` using `active_spins=find(opspec)`.
- Lines 57: computes `[rows,cols,vals]` using `[rows,cols,vals]=find(A)`.
- Lines 63: computes `source` using `source=cell(1,numel(active_spins))`.
- Lines 66: computes `destin` using `destin=cell(1,numel(active_spins))`.
- Lines 69: computes `struct` using `struct=cell(1,numel(active_spins))`.
- Lines 75: computes `mult` using `mult=spin_system.comp.mults(active_spins(n))`.
- Lines 76: computes `table_idx` using `table_idx=opspec(active_spins(n))+1`.
- Lines 84: computes `pt` using `pt=squeeze(spin_system.bas.lpst{mult}(table_idx,:,:))`.
- Lines 99: computes `[destin{n},source{n},struct{n}]` using `[destin{n},source{n},struct{n}]=find(pt)`.
- Lines 102: computes `source{n}` using `source{n}=source{n}-1; destin{n}=destin{n}-1`.
- Lines 107: computes `from` using `from=source{1}; to=destin{1}; coeff=struct{1}`.
- Lines 111-112: computes `to` using `to=[kron(to,ones(size(destin{n},1),1)) kron(ones(size(to,1),1),destin{n})]`.
- Lines 113: computes `coeff` using `coeff=kron(coeff,struct{n})`.
- Lines 117: computes `basis_cols` using `basis_cols=spin_system.bas.basis(:,active_spins)`.
- Lines 122: computes `from(kill_mask,:)` using `from(kill_mask,:)=[]; to(kill_mask,:)=[]; coeff(kill_mask,:)=[]`.

### Local helper functions

- Line 188: `grumble()` — `function grumble(spin_system,opspec)`.
  - Representative operation: `if ~isfield(spin_system,'bas')`.
  - Representative operation: `error('basis set information is missing, run basis() before calling this function.')`.

## Parameters / inputs

- opspec -Spinach operator specification described in Sections 2.1
- and 3.3 of the following paper:
- side -'left' or 'right' causes the function to return a product
- superoperator corresponding to a product from that side;
- 'comm' or 'acomm' results in commutation and anticommuta-
- tion superoperator respectively.

## Outputs

- A -a three-column array of row indices (first column),
- column indices (second column) and values (third column).
- Note: this is a very general function to which direct calls are not
- usually required -please use the (much friendlier) operator()
- function.
- Note: the superoperator is returned in XYZ sparse format, which is
- different from Matlab's CSC format.

## Implementation structure

- Sided product superoperator in the spherical tensor basis set. Returns
- superoperators corresponding to right or left multiplication of a den-
- sity matrix by a user-specified operator. Syntax:
- A=superop(spin_system,opspec,side)
- opspec -Spinach operator specification described in Sections 2.1
- and 3.3 of the following paper:
- side -'left' or 'right' causes the function to return a product
- superoperator corresponding to a product from that side;
- 'comm' or 'acomm' results in commutation and anticommuta-
- tion superoperator respectively.
- A -a three-column array of row indices (first column),
- column indices (second column) and values (third column).

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `strcmp()`, `grumble()`, `unit_oper()`, `active_spins()`, `opspec()`, `squeeze()`, `ismember()`, `from()`, `coeff()`, `true()`, `and()`, `basis_cols()`, `source_subsp()`, `destin_subsp()`, `isequal()`, `source_subsp_idx()`.
