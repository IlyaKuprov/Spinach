# kernel/optimcon/alpha_conds.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/optimcon/alpha_conds.m`
- Signature: `test=alpha_conds(test_type,alpha,fx_0,fx_1,gfx_0,gfx_1,dir,spin_system)`
- Total lines: 125

## Purpose

Applies one of the line search acceptance tests used by the brac- keting and sectioning routines in constrained optimisation and returns true when the chosen condition is satisfied. Syntax: test=alpha_conds(test_type,alpha,fx_0,fx_1,... gfx_0,gfx_1,dir,spin_system)

## Physical / mathematical content

- Optimal-control core routines. These files implement GRAPE-style objective evaluation, quasi-Newton search, line search, regularisation, distortion models, and waveform parameterisations.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 44-45: Check consistency; implemented by `grumble(test_type,alpha,fx_0,fx_1,gfx_0,gfx_1,dir)`.
- Lines 47-48: Test selection; implemented by `if test_type==0`.
- Lines 50-51: Monotonic increase; implemented by `test=(fx_1>fx_0)`.
- Lines 55-56: Armijo sufficient increase condition; implemented by `test=(fx_1>=fx_0+spin_system.control.ls_c1*alpha*(gfx_0'*dir))`.
- Lines 60-61: Strong Wolfe curvature condition; implemented by `test=(abs(gfx_1'*dir)<=spin_system.control.ls_c2*abs(gfx_0'*dir))`.
- Lines 65-67: Positive derivative in descent direction; implemented by `test=(gfx_1'*dir>0)`.

### Control flow inferred from the code

- Line 48: conditional branch on `test_type==0`.

### Key state/data transformations

- Lines 51: computes `test` using `test=(fx_1>fx_0)`.

### Local helper functions

- Line 74: `grumble()` — `function grumble(test_type,alpha,fx_0,fx_1,gfx_0,gfx_1,dir)`.
  - Representative operation: `if (~isnumeric(test_type))||(~isreal(test_type))||(~isscalar(test_type))|| (~ismember(test_type,[0 1 2 3]))`.
  - Representative operation: `(~ismember(test_type,[0 1 2 3]))`.

## Parameters / inputs

- test_type -condition selector:
- 0 for monotonic increase test
- 1 for Armijo sufficient increase test
- 2 for strong Wolfe curvature test
- 3 for ascent direction test
- alpha -trial step length
- fx_0 -objective value at the initial point
- fx_1 -objective value at the trial point
- gfx_0 -gradient at the initial point
- gfx_1 -gradient at the trial point
- dir -search direction vector
- spin_system -Spinach data structure with line
- search settings in control

## Outputs

- test -logical true if the selected
- condition is satisfied

## Implementation structure

- Applies one of the line search acceptance tests used by the brac-
- keting and sectioning routines in constrained optimisation and
- returns true when the chosen condition is satisfied. Syntax:
- test=alpha_conds(test_type,alpha,fx_0,fx_1,...
- gfx_0,gfx_1,dir,spin_system)
- test_type -condition selector:
- 0 for monotonic increase test
- 1 for Armijo sufficient increase test
- 2 for strong Wolfe curvature test
- 3 for ascent direction test
- alpha -trial step length
- fx_0 -objective value at the initial point

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `isscalar()`, `ismember()`, `iscolumn()`, `isequal()`.
