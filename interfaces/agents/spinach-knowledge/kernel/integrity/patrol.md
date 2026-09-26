# kernel/integrity/patrol.m

- Signature: `patrol(test_subject)`

## Purpose

This function runs contiouously on one of our servers, its purpose is to catch any unintended consequences before they propagate too far down the development chain. Syntax: patrol(test_subject)

## Physical / mathematical content

- Integrity-control utilities. These files check distribution state, path collisions, style conformance, sniffer databases, and other safeguards that protect Spinach reproducibility.

## Numerical / algorithmic content

## Parameters / inputs

- test_subject -a character string; if it occurs
- anywhere within the example file
- path, that file is included into
- the patrol run

## Outputs

- whatever the individual examples return

## Implementation structure

- This function runs contiouously on one of our servers, its
- purpose is to catch any unintended consequences before they
- propagate too far down the development chain. Syntax:
- patrol(test_subject)
- test_subject -a character string; if it occurs
- anywhere within the example file
- path, that file is included into
- the patrol run
- whatever the individual examples return
- Set default and check consistency
- List exceptions
- Shuffle the RNG
