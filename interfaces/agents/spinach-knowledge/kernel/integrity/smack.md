# kernel/integrity/smack.m

[MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/integrity/smack.m) · [Spin Dynamics Wiki](https://spindynamics.org/wiki/index.php?title=smack.m)

`smack()` is a MATLAB environment-recovery utility for problems involving MATLAB Distributed Computing Server (MDCS), not a spin-dynamics calculation. Its source comment says to use it from the command line; the function does not enforce that restriction.

The function deletes the current parallel pool, deletes jobs belonging to the `Processes` cluster, closes all open MATLAB file handles, clears the workspace, and resets each device counted by `gpuDeviceCount`. These are broad session side effects, not a selective cleanup of a particular job or GPU.

## Inputs and outputs

No inputs or return value. There is no source-level validation or recovery guard around the cleanup operations. Use only when those whole-session effects are intended.
