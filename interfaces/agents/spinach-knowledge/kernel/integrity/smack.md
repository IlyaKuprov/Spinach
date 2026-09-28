# kernel/integrity/smack.m

- Signature: `smack()`

## Purpose

Resets Matlab state after problems with MDCS (Matlab Distributed Computing Server). This function is intended for use from the command line only.

## Physical / mathematical content

This is a Matlab environment-recovery utility; it does not perform a physical calculation.

## Numerical / algorithmic content

No numerical calculation is performed. The function deletes the current parallel pool and jobs on the `Processes` cluster, closes open file handles, clears the workspace, and resets available GPUs.

## Parameters / inputs

None.

## Outputs

No return value. Matlab state is cleared and any available GPU devices are reset.

## Implementation structure

The function deletes `gcp('nocreate')`, deletes jobs from `parcluster('Processes')`, calls `fclose('all')` and `clear('all')`, then resets each device returned by `gpuDeviceCount`.
