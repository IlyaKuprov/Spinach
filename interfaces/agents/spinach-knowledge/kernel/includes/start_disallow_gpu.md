# kernel/includes/start_disallow_gpu.m

- Signature: `(script file)`

## Purpose

Temporarily disable GPU use when it is enabled in `spin_system.sys.enable`, preserving the previous state for `end_disallow_gpu`.

## Behaviour

The include records whether `'gpu'` is present in `spin_system.sys.enable` as `user_wanted_gpu` and always reports a programmer-requested GPU-disallow warning. If GPU use was enabled, it removes `'gpu'` from the enable list.

## Source documentation

https://spindynamics.org/wiki/index.php?title=start_disallow_gpu.m
