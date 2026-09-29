# kernel/integrity/existentials.m

Source: [MATLAB implementation](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/integrity/existentials.m) · [Spin Dynamics Wiki](https://spindynamics.org/wiki/index.php?title=existentials.m).

- Signature: `existentials()`

## Purpose

A startup integrity check for the MATLAB installation and Spinach path. It returns without checking when called on a parallel worker; on the client it reports that startup checks are running, then verifies environment prerequisites and the visibility of Spinach files.

## Integrity lifecycle

The client-side checks proceed in source order: MATLAB must be R2026a or newer; the MATLAB installation must contain Parallel Computing, Deep Learning, Reinforcement Learning, Optimisation, Statistics and Machine Learning, and Mapping toolboxes; then the routine recursively enumerates Spinach `.m` files under `kernel`, `interfaces`, `experiments`, and `etc`.

For each enumerated file it constructs the expected Spinach pathname and asks MATLAB `which` for the resolved file. A different resolved pathname is reported as a same-name collision, except for overloads; an empty result is reported as a Spinach file missing from the MATLAB path. Those path failures print diagnostic paths and stop with `startup checks not passed`. On Windows, for a drive-letter path, a non-NTFS filesystem produces a warning about reliable file locking and performance.

## Inputs, outputs, and units

There are no function arguments or returned values. This is an environment/path check, not a numerical or physical calculation; equations, normalisation, matrix shape, and Hz-versus-angular-frequency units do not apply.

## Source guards

The worker early return is `if isworkernode, return; end`. The routine is not parameterised by a user-supplied input domain. Missing release/toolbox prerequisites fail immediately with explicit errors before the recursive path audit.

Related source-backed checks: [`exorcise.m`](./exorcise.md) audits source conventions; [`patrol.m`](./patrol.md) selects and runs examples.
