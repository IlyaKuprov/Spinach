# kernel/integrity/sniff.m

- Signature: `sniff(action)`

## Purpose

Kernel integrity control. Checks Spinach distribution .m files for any modifications that the user did since downloading Spi- nach. The function prints the list of files that have changed in any way since the internal database has been rearmed. The purpose is to catch local modifications that the user may have made and forgotten about, that are causing some unintend- ed consequences elsewhere in Spinach. Syntax: snif

## Physical / mathematical content

- Integrity-control utilities. These files check distribution state, path collisions, style conformance, sniffer databases, and other safeguards that protect Spinach reproducibility.

## Numerical / algorithmic content

## Parameters / inputs

- action -'none' prints the names of fishy files to
- the console, 'open' opens them

## Implementation structure

- Kernel integrity control. Checks Spinach distribution .m files
- for any modifications that the user did since downloading Spi-
- nach. The function prints the list of files that have changed
- in any way since the internal database has been rearmed.
- The purpose is to catch local modifications that the user may
- have made and forgotten about, that are causing some unintend-
- ed consequences elsewhere in Spinach. Syntax:
- sniff(action)
- action -'none' prints the names of fishy files to
- the console, 'open' opens them
- Default is to take no action
- Check consistency
