# kernel/includes/end_disallow_gpu.m

- Signature: script include; reads variables in the shared script workspace

## Purpose and behaviour

This is the closing include for the GPU-disallow sequence. It requires the variable `user_wanted_gpu`; if that variable does not exist, it errors with a message requiring a preceding `start_disallow_gpu` command. When `user_wanted_gpu` is true, it appends `'gpu'` to `spin_system.sys.enable`. When false, it makes no change. The source does not otherwise restore or validate the enable list, nor does it check for an existing duplicate before appending.

## Source comment attribution

The source quotes Юлий Ким, “Истерическая перестроечная” (1988).

## Source links

- [MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/includes/end_disallow_gpu.m)
- [Spinach Wiki](https://spindynamics.org/wiki/index.php?title=end_disallow_gpu.m)
