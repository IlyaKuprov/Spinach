# kernel/includes/autoexec.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/includes/autoexec.m`
- Signature: `(script file)`
- Total lines: 71

## Purpose

This include is executed at the start of create.m, it over- rides all user input. A good use case is forcing polyadic or GPU arithmetic, or some other specific hardware or soft- ware configuration.

## Physical / mathematical content

- Include scripts and shared setup fragments. These files implement tightly scoped runtime setup, parallel profiling, resource guards, or shared kernels included by other Spinach routines.

## Numerical / algorithmic content


## Code-derived implementation details

### Comment-guided execution stages

- Lines 1-11: This include is executed at the start of create.m, it over- rides all user input. A good use case is forcing polyadic or GPU arithmetic, or some other specific hardware or soft- ware configuration.; implemented by `warning('off','parallel:gpu:device:DeviceDeprecated')`.
- Lines 10-11: Kill the pointless GPU deprecation warning; implemented by `warning('off','parallel:gpu:device:DeviceDeprecated')`.
- Lines 13-14: Kill stupid ass figure defaults in R2025a and later; implemented by `set(groot,'defaultFigurePosition',[680 458 560 420])`.
- Lines 19-20: Do not override user spec; implemented by `if ~isfield(sys,'parallel')`.
- Lines 22-23: IK group system settings; implemented by `switch getenv('COMPUTERNAME')`.
- Lines 27-29: Be careful with GPUs; implemented by `if isfield(sys,'enable')&& ismember('gpu',sys.enable)`.
- Lines 31-32: 4 workers per GPU are safe; implemented by `sys.parallel={'processes',32}`.
- Lines 42-43: 4 workers per GPU are safe; implemented by `sys.parallel={'processes',12}`.
- Lines 49-51: Do nothing; implemented by `end`.

### Control flow inferred from the code

- Line 20: conditional branch on `~isfield(sys,'parallel')`.
- Line 23: dispatches on `getenv('COMPUTERNAME')`; cases `'ALAUNDO'`, `'TALOS'`.
- Line 28: conditional branch on `isfield(sys,'enable')&&`.
- Line 39: conditional branch on `isfield(sys,'enable')&&`.

### Key state/data transformations

- Lines 32: computes `sys.parallel` using `sys.parallel={'processes',32}`.

## Implementation structure

- This include is executed at the start of create.m, it over-
- rides all user input. A good use case is forcing polyadic
- or GPU arithmetic, or some other specific hardware or soft-
- ware configuration.
- Kill the pointless GPU deprecation warning
- Kill stupid ass figure defaults in R2025a and later
- Do not override user spec
- IK group system settings
- Be careful with GPUs
- 4 workers per GPU are safe
- Do nothing
- This relocates the scratch folder

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `set()`, `isfield()`, `getenv()`, `ismember()`.
