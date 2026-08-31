# kernel/residual.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/residual.m`
- Signature: `spin_system=residual(spin_system)`
- Total lines: 106

## Purpose

Sets up interaction tensors under partial ordering in a liquid crystal with the user-supplied order matrix. All adjustable pa- rameters are set during the call to create.m function. Syntax: spin_system=residual(spin_system)

## Physical / mathematical content

- This file belongs to the `kernel` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 38-39: Check consistency; implemented by `grumble(spin_system)`.
- Lines 41-42: Loop over chemical subsystems; implemented by `for s=1:numel(spin_system.chem.parts)`.
- Lines 44-45: Process Zeeman interactions; implemented by `for n=spin_system.chem.parts{s}`.
- Lines 47-48: Process non-empty interaction tensors; implemented by `if ~isempty(spin_system.inter.zeeman.matrix{n})`.
- Lines 50-51: Obtain isotropic part; implemented by `iso=eye(3)*trace(spin_system.inter.zeeman.matrix{n})/3`.
- Lines 53-54: Calculate residual order; implemented by `extra_zz=trace(spin_system.inter.order_matrix{s}*(spin_system.inter.zeeman.matrix{n}-iso))`.
- Lines 56-57: Update Zeeman tensor; implemented by `spin_system.inter.zeeman.matrix{n}=iso+diag([-extra_zz/3 -extra_zz/3 2*extra_zz/3])`.
- Lines 63-64: Process spin-spin couplings; implemented by `for n=spin_system.chem.parts{s}`.
- Lines 67-68: Process non-empty interaction tensors; implemented by `if ~isempty(spin_system.inter.coupling.matrix{n,k})`.
- Lines 70-71: Obtain isotropic part; implemented by `iso=trace(spin_system.inter.coupling.matrix{n,k})*eye(3)/3`.
- Lines 73-74: Calculate residual order; implemented by `extra_zz=trace(spin_system.inter.order_matrix{s}*(spin_system.inter.coupling.matrix{n,k}-iso))`.
- Lines 76-77: Update coupling tensor; implemented by `spin_system.inter.coupling.matrix{n,k}=iso+diag([-extra_zz/3 -extra_zz/3 2*extra_zz/3])`.
- Lines 81-82: Drop the couplings that become insignificant after averaging; implemented by `if norm(spin_system.inter.coupling.matrix{n,k},2)<2*pi*spin_system.tols.inter_cutoff`.
- Lines 91-92: Report back to the user; implemented by `report(spin_system,'interaction tensors have been replaced by their weak order residuals.')`.

### Control flow inferred from the code

- Line 42: `for` loop over `s=1:numel(spin_system.chem.parts)`.
- Line 45: `for` loop over `n=spin_system.chem.parts{s}`.
- Line 48: conditional branch on `~isempty(spin_system.inter.zeeman.matrix{n})`.
- Line 64: `for` loop over `n=spin_system.chem.parts{s}`.
- Line 65: `for` loop over `k=spin_system.chem.parts{s}`.
- Line 68: conditional branch on `~isempty(spin_system.inter.coupling.matrix{n,k})`.
- Line 82: conditional branch on `norm(spin_system.inter.coupling.matrix{n,k},2)<2*pi*spin_system.tols.inter_cutoff`.

### Key state/data transformations

- Lines 51: computes `iso` using `iso=eye(3)*trace(spin_system.inter.zeeman.matrix{n})/3`.
- Lines 54: computes `extra_zz` using `extra_zz=trace(spin_system.inter.order_matrix{s}*(spin_system.inter.zeeman.matrix{n}-iso))`.
- Lines 57: computes `spin_system.inter.zeeman.matrix{n}` using `spin_system.inter.zeeman.matrix{n}=iso+diag([-extra_zz/3 -extra_zz/3 2*extra_zz/3])`.
- Lines 77: computes `spin_system.inter.coupling.matrix{n,k}` using `spin_system.inter.coupling.matrix{n,k}=iso+diag([-extra_zz/3 -extra_zz/3 2*extra_zz/3])`.

### Local helper functions

- Line 97: `grumble()` — `function grumble(spin_system)`. If I'd observed all the rules, I'd never have got anywhere. Marilyn Monroe
  - Representative operation: `if (~isfield(spin_system.inter,'order_matrix'))||isempty(spin_system.inter.order_matrix)`.
  - Representative operation: `error('order matrix infomation is missing from the spin_system structure.')`.

## Parameters / inputs

- spin_system -the output of create.m containing
- spin system and interaction infor-
- mation, which must include the or-
- der matrix

## Outputs

- spin_system -the same object with anisotropic
- parts of all interaction tensors
- replaced with their partial order
- residuals.
- Note: this function is only applicable to weak residual order
- in high-field NMR spectroscopy.
- Note: the function overwrites the interaction tensors supplied
- by the user. Relaxation superoperator, if required, must
- be computed before this function is called.
- Note: this function is invoked automatically by liquid.m con-
- text when when parameters.needs cell array contains 'rdc'.

## Implementation structure

- Sets up interaction tensors under partial ordering in a liquid
- crystal with the user-supplied order matrix. All adjustable pa-
- rameters are set during the call to create.m function. Syntax:
- spin_system=residual(spin_system)
- spin_system - the output of create.m containing
- spin system and interaction infor-
- mation, which must include the or-
- der matrix
- spin_system - the same object with anisotropic
- parts of all interaction tensors
- replaced with their partial order
- residuals.

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `report()`, `isfield()`.
