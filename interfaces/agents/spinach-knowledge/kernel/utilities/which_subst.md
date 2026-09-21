# kernel/utilities/which_subst.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/utilities/which_subst.m`
- Signature: `subst=which_subst(spin_system,spins)`
- Total lines: 62

## Purpose

Finds out which substance hosts the specified spins; throws an error if there is more than one. Syntax: subst=which_subst(spin_system,spins)

## Physical / mathematical content

- General mathematical and infrastructure utilities. This area contains finite differences, perturbation theory, graph algorithms, spectral densities, tensor algebra, hash/report helpers, and other reusable numerical components.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Parameters / inputs

- spins -a list of positive integers

## Outputs

- subst -a positive integer

## Implementation structure

- Finds out which substance hosts the specified spins;
- throws an error if there is more than one. Syntax:
- subst=which_subst(spin_system,spins)
- spins -a list of positive integers
- subst -a positive integer
- Check consistency
- Find the substances hosting specified spins
- Only one substance is permitted
- Get substance number
- Confirm that all spins are in the same substance
- Consistency enforcement
- I swear to you that to think too much is

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `cellfun()`, `any()`, `ismember()`, `nnz()`, `all()`, `isvector()`.
