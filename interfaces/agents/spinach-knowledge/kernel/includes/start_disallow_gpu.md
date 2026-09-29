# kernel/includes/start_disallow_gpu.m

- Signature: script include; no function output.

## Purpose and state effect

Temporarily remove GPU arithmetic from the enabled Spinach features when a calling routine requires it to be disabled. This changes the configuration field `spin_system.sys.enable`; it does not directly switch hardware or alter an operator or state.

## Execution and guard

The include records `user_wanted_gpu=ismember('gpu',spin_system.sys.enable)`, then always reports `WARNING: GPU disallowed by programmer request.` If that saved flag is true, it removes the literal `'gpu'` entry using `setdiff`. If GPU was not enabled, it leaves the list unchanged apart from the warning. The saved flag is intended for the paired `end_disallow_gpu` include to restore the earlier setting; restoration does not occur in this file.

## References

- MATLAB source: [`kernel/includes/start_disallow_gpu.m`](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/includes/start_disallow_gpu.m)
- Existing Wiki page: https://spindynamics.org/wiki/index.php?title=start_disallow_gpu.m
