# kernel/includes/start_disallow_gpu.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/includes/start_disallow_gpu.m`
- Signature: `(script file)`
- Total lines: 25

## Purpose

Forces GPU arithmetic to be turned off even if the user had requested it in sys.enable setting; restore previous state using end_disallow_gpu command.

## Physical / mathematical content

- Include scripts and shared setup fragments. These files implement tightly scoped runtime setup, parallel profiling, resource guards, or shared kernels included by other Spinach routines.

## Numerical / algorithmic content


## Code-derived implementation details

### Comment-guided execution stages

- Lines 1-10: Forces GPU arithmetic to be turned off even if the user had requested it in sys.enable setting; restore previous state using end_disallow_gpu command.; implemented by `user_wanted_gpu=ismember('gpu',spin_system.sys.enable)`.
- Lines 9-10: Check if GPU is currently enabled; implemented by `user_wanted_gpu=ismember('gpu',spin_system.sys.enable)`.
- Lines 13-14: Disable GPU if it had been enabled; implemented by `if user_wanted_gpu`.

### Control flow inferred from the code

- Line 14: conditional branch on `user_wanted_gpu`.

### Key state/data transformations

- Lines 10: computes `user_wanted_gpu` using `user_wanted_gpu=ismember('gpu',spin_system.sys.enable)`.
- Lines 15: computes `spin_system.sys.enable` using `spin_system.sys.enable=setdiff(spin_system.sys.enable,{'gpu'})`.

## Implementation structure

- Forces GPU arithmetic to be turned off even if the user had
- requested it in sys.enable setting; restore previous state
- using end_disallow_gpu command.
- Check if GPU is currently enabled
- Disable GPU if it had been enabled
- Once at MIT, I saw a frat bro accidentally smear a line all over the
- table. Somehow he still managed to snort the whole thing. He looked
- me square in the eyes and said "same high bro, Stokes theorem".
- Internet folklore
- #NHEAD #NGRUM

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `ismember()`, `report()`, `setdiff()`.
