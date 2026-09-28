# kernel/integrity/rearm.m

- Signature: `rearm()`

## Purpose

Rebuilds the sniffer database from the current contents of Spinach's `.m` files. The `sniff` function uses this baseline to identify files changed since it was created.

## Physical / mathematical content

This is an integrity utility; it does not model a physical system.

## Numerical / algorithmic content

For each included file, it hashes the filename together with a hash of the file contents. The resulting list is saved as `smells.mat`, replacing the previous database.

## Parameters / inputs

None.

## Outputs

No return value. The function overwrites `smells.mat` and displays `rearm: sniffer rearmed.`

## Implementation structure

It scans `.m` files under `kernel`, `interfaces`, `experiments`, and `etc`, applies the exception list, collects the hashes, deletes the existing `smells.mat`, and saves the new `smells` variable.
