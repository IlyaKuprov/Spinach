# kernel/includes/autoexec.m

- Signature: `(script file)`

## Purpose

This include runs at the start of `create.m`. It disables the GPU device deprecation warning and sets default figure position, window style, menu bar, and toolbar. It does not replace a user-specified `sys.parallel` value.

## Behaviour

When `sys.parallel` is absent, GPU use is enabled, and the host is `ALAUNDO` or `TALOS`, the include sets process parallelism to 32 or 12, respectively. The scratch-folder assignment in the source is commented out and has no effect.

## Source documentation

https://spindynamics.org/wiki/index.php?title=autoexec.m
