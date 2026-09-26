# kernel/residual.m

- Signature: `spin_system=residual(spin_system)`

## Purpose

Sets up interaction tensors under partial ordering in a liquid crystal with the user-supplied order matrix. All adjustable pa- rameters are set during the call to create.m function. Syntax: spin_system=residual(spin_system)

## Physical / mathematical content

## Numerical / algorithmic content

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
