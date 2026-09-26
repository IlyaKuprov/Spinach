# kernel/includes/start_disallow_gpu.m

- Signature: `(script file)`

## Purpose

Forces GPU arithmetic to be turned off even if the user had requested it in sys.enable setting; restore previous state using end_disallow_gpu command.

## Physical / mathematical content

- Include scripts and shared setup fragments. These files implement tightly scoped runtime setup, parallel profiling, resource guards, or shared kernels included by other Spinach routines.

## Numerical / algorithmic content

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
