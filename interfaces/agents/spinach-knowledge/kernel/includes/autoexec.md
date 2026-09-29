# kernel/includes/autoexec.m

- Signature: script include, executed at the start of `create.m`

## Purpose and unconditional setup

The include turns off the `parallel:gpu:device:DeviceDeprecated` warning and sets root-graphics defaults: figure position `[680 458 560 420]`, normal window style, figure menu bar, and figure toolbar. The source comment frames the figure settings as overrides for MATLAB R2025a and later defaults. These statements run independently of the host and parallel settings.

## Host and parallel-setting guards

Host-specific setup runs only if `sys.parallel` does not already exist. The host switch uses the exact value of `getenv('COMPUTERNAME')`:

- On `ALAUNDO` (source comment: 128 Intel cores, 4 TB RAM, 8 H200 GPUs), it sets `sys.parallel={'processes',32}` only when `sys.enable` exists and contains `'gpu'`. The source comment describes four workers per GPU as safe.
- On `TALOS` (source comment: 56 Intel cores, 1 TB RAM, 3 A800 GPUs), it sets `sys.parallel={'processes',12}` only when `sys.enable` exists and contains `'gpu'`, again four workers per GPU by the source comment.
- On any other host, or on either named host without that GPU-enabled condition, the switch does not assign `sys.parallel`.

Thus an existing user value is preserved, and a missing value is not automatically filled on CPU-only configurations. The apparent scratch relocation is commented-out example text and has no runtime effect.

## Source links

- [MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/includes/autoexec.m)
- [Spinach Wiki](https://spindynamics.org/wiki/index.php?title=autoexec.m)
