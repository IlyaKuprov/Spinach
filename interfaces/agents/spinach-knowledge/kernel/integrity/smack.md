# kernel/integrity/smack.m

- Signature: `smack()`

## Purpose

Gives Matlab a good smack every time MDCS gets its kni- ckers in a twist. Syntax: smack() This function shuts down the parallel pool, clears the workspace, clears the GPUs, and makes sure there are no crashed MDCS jobs left over. This function should only be used from the command line.

## Physical / mathematical content

- Integrity-control utilities. These files check distribution state, path collisions, style conformance, sniffer databases, and other safeguards that protect Spinach reproducibility.

## Numerical / algorithmic content

## Implementation structure

- Gives Matlab a good smack every time MDCS gets its kni-
- ckers in a twist. Syntax:
- smack()
- This function shuts down the parallel pool, clears the
- workspace, clears the GPUs, and makes sure there are no
- crashed MDCS jobs left over. This function should only
- be used from the command line.
- Kill the parallel pool
- Clear out crashed jobs
- Close all handles
- Clear the workspace
- Reset all GPUs
