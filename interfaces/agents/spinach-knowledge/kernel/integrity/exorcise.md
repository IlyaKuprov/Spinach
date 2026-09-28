# kernel/integrity/exorcise.m

- Signature: `exorcise(mode)`

## Purpose

Scans Spinach source files for violations of the repository's house style and opens the first offending file before reporting an error.

## Physical / mathematical content

This is a source-integrity utility; it does not perform a physical or numerical calculation.

## Numerical / algorithmic content

The scan checks file formatting and documentation, MATLAB syntax, and several coding conventions. In `online` mode it also checks that each documented Wiki page is available; `offline` mode skips that network check.

## Parameters / inputs

- `mode` — `'online'` checks the corresponding documentation Wiki page; `'offline'` skips the Wiki check.

## Outputs

No return value. The function reports scan progress and success; on the first detected violation it opens the file in the editor and raises an error.

## Implementation structure

It visits `.m` files under `kernel`, `interfaces`, `experiments`, and `etc` in randomized order, excluding the `jsonlab-1.5` foreign-package directory. Checks include required headers and `grumble` validation, whitespace and tab rules, MATLAB's `checkcode`, explicit norm types, portable path separators, a top-level `otherwise` in each `switch`, and use of `report` rather than `disp` when `spin_system` is available.
