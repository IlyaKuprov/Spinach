# kernel/integrity/existentials.m

- Signature: `existentials()`

## Purpose

Checks the Matlab environment and Spinach path during startup. It returns immediately on parallel workers. On the client, it checks the Matlab release and required toolboxes, then detects function-name collisions and Spinach files that are not visible on the Matlab path.

## Physical / mathematical content

This is an environment-integrity check; it does not model a physical system.

## Numerical / algorithmic content

No numerical calculation is performed. The routine compares each discovered Spinach file with the location returned by Matlab's `which`.

## Parameters / inputs

None.

## Outputs

No return value. It displays startup-check progress and raises an error if a prerequisite, path entry, or collision check fails.

## Implementation structure

The routine requires Matlab R2026a or later and the Parallel Computing, Deep Learning, Reinforcement Learning, Optimisation, Statistics and Machine Learning, and Mapping toolboxes. It scans `.m` files under `kernel`, `interfaces`, `experiments`, and `etc`. A same-named file outside Spinach is reported as a collision, except for overloads; a Spinach file that `which` cannot find is reported as a path setup problem.
