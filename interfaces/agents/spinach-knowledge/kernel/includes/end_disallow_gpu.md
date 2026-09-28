# kernel/includes/end_disallow_gpu.m

- Signature: `(script file)`

## Purpose

Restore GPU availability after `start_disallow_gpu.m` when GPU use had been enabled before that include ran.

## Behaviour

The include errors if `user_wanted_gpu` is absent. If it is true, the include appends `'gpu'` to `spin_system.sys.enable`; otherwise it leaves the GPU setting unchanged.

## Source comment attribution

Юлий Ким, "Истерическая перестроечная", 1988.

## Source documentation

https://spindynamics.org/wiki/index.php?title=end_disallow_gpu.m
